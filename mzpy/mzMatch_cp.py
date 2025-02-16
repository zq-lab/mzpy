import cupy as cp
import numpy as np
import pandas as pd


class Match():
    def __init__(  
                self,  
                ref_mz,  
                ref_ms,  
                que_mz=None,  
                que_ms=None,  
                sort_ms=False 
    ):  
        """  
        优化版 GPUmatch:  
          1) _cosine_ 使用 searchsorted 匹配，而非 Python BFS 循环  
          2) 仅对满足母离子公差的谱对调用 _cosine_，减少不必要计算  
          3) 尽量将大规模操作放在 GPU 上做向量化  
        """  

        # 如果 que_mz, que_ms 为空，则认为自对比  
        self.is_symmetric = (que_mz is None or len(que_mz) == 0) and (que_ms is None or len(que_ms) == 0)  

        # 转到 GPU (float32)  
        self.ref_mz = cp.asarray(ref_mz, dtype=cp.float32)
        self.ref_ms = []  
        for item in ref_ms:  
            arr = cp.asarray(item, dtype=cp.float32)  
            if sort_ms:  
                idx = cp.argsort(arr[:,0])  
                arr = arr[idx]  
            self.ref_ms.append(arr)  
            
        if self.is_symmetric:
            self.que_mz = None
            self.que_ms = None
        else:
            self.que_mz = cp.asarray(que_mz, dtype=cp.float32)  
            self.que_ms = []  
            for item in que_ms:  
                arr = cp.asarray(item, dtype=cp.float32)  
                if sort_ms:  
                    idx = cp.argsort(arr[:,0])  
                    arr = arr[idx]  
                self.que_ms.append(arr)        
        
    def _cosine_(self, a, b, tol=0.005):
        """
        在 GPU 上计算两条已在 GPU 上的 msms 谱 a, b 的“匹配率修正”余弦相似度。
        a, b: shape (N,2), (M,2)，分别表示 [mz, intensity]，且均为 cupy.ndarray。
        tol: 峰匹配的公差 (绝对差值)
    
        算法:
        1) 先按 mz 升序排序两条谱。
        2) 逐峰比对，当 |mz_a - mz_b| <= tol 视为匹配并记录。
        3) 计算普通余弦相似度 cos_val = dot_val / (norm_a * norm_b)。
        4) 匹配率 = matched_count / (n_a + n_b)。
        5) 最终返回 cos_val * 匹配率。
        """    
        # 按 mz 升序排序
        sort_idx_a = cp.argsort(a[:, 0])
        sort_idx_b = cp.argsort(b[:, 0])
        a_sorted = a[sort_idx_a]
        b_sorted = b[sort_idx_b]
    
        n_a = a_sorted.shape[0]
        n_b = b_sorted.shape[0]
    
        # 准备并集数组(最大长度=n_a+n_b)
        union_a = cp.empty(n_a + n_b, dtype=a.dtype)
        union_b = cp.empty(n_a + n_b, dtype=b.dtype)
    
        i = 0
        j = 0
        idx = 0
        matched_count = 0  # 用于统计峰匹配数量
    
        while i < n_a and j < n_b:
            mz_a, int_a = a_sorted[i]
            mz_b, int_b = b_sorted[j]
            if cp.abs(mz_a - mz_b) <= tol:
                union_a[idx] = int_a
                union_b[idx] = int_b
                i += 1
                j += 1
                matched_count += 1  # 匹配峰数量+1
            elif mz_a < mz_b:
                union_a[idx] = int_a
                union_b[idx] = 0.0
                i += 1
            else:
                union_a[idx] = 0.0
                union_b[idx] = int_b
                j += 1
            idx += 1
    
        # 补齐剩余
        while i < n_a:
            union_a[idx] = a_sorted[i, 1]
            union_b[idx] = 0.0
            i += 1
            idx += 1
    
        while j < n_b:
            union_a[idx] = 0.0
            union_b[idx] = b_sorted[j, 1]
            j += 1
            idx += 1
    
        # 截断到实际并集长度
        union_a = union_a[:idx]
        union_b = union_b[:idx]
    
        # 计算余弦
        norm_a = cp.sqrt(cp.sum(union_a * union_a))
        norm_b = cp.sqrt(cp.sum(union_b * union_b))
        if norm_a == 0 or norm_b == 0:
            return 0.0
    
        dot_val = cp.sum(union_a * union_b)
        cos_val = dot_val / (norm_a * norm_b)
    
        # 匹配率修正
        match_ratio = matched_count / idx if idx > 0 else 0.0
        return cos_val * match_ratio
    

    def cosine_mx(self, tol=(0.003, 0.005)):
        """
        计算两批母离子 ref_mz, query_mz 的匹配后 msms 余弦相似度矩阵 (匹配率修正版本)，GPU加速。
    
        tol=(match_tol, msms_tol):
           tol[0] -> 母离子 mz 匹配公差
           tol[1] -> msms 峰比对时的 mz 公差
    
        返回
        ----------
        sim_matrix : cupy.ndarray, shape (R, Q)
            sim_matrix[i, j] = 当 |ref_mz[i] - query_mz[j]| <= tol[0] 时，计算匹配率修正的余弦相似度；
                               否则为 0.0。
        """
        if self.is_symmetric:
            que_mz = self.ref_mz #这是临时赋值，不要去更改self.que_mz and self.que_ms
            que_ms = self.ref_ms
        else:
            que_mz = self.que_mz
            que_ms = self.que_ms

        R = self.ref_mz.shape[0]  
        Q = que_mz.shape[0]  
        similarity = cp.full((R, Q), -1.0, dtype=cp.float32) 
        # 广播比对母离子  
        diff_mz = cp.abs(self.ref_mz[:, None] - que_mz[None, :])  
        matched_mask = diff_mz <= tol[0]
        matched_i, matched_j = cp.where(matched_mask)  

        if self.is_symmetric:  
            # 只保留上三角(含对角线)  
            upper_mask = matched_j >= matched_i  
            matched_i = matched_i[upper_mask]  
            matched_j = matched_j[upper_mask]  

            # 对角线初始化为 1.0  
            diag_idx = cp.arange(R)  
            similarity[diag_idx, diag_idx] = 1.0  

            # 遍历 matched_i, matched_j 完成上三角  
            for k in range(matched_i.shape[0]):  
                i = matched_i[k].item()  
                j = matched_j[k].item()  
                if i == j:  
                    continue  # 对角线已是 1.0, 需要可改成真实计算则移除此行  
                val = self._cosine_(self.ref_ms[i], que_ms[j], tol[1])  
                similarity[i, j] = val  
                similarity[j, i] = val  # 镜像  
        else:  
            for k in range(matched_i.shape[0]):  
                i = matched_i[k].item()  
                j = matched_j[k].item()  
                val = self._cosine_(self.ref_ms[i], que_ms[j], tol[1])  
                similarity[i, j] = val  
    
        return similarity.get() #从GPU返回

    @classmethod
    def cosine_df_chunk(cls, df, num_chunks, csv_fname, mz_on='precursormz', ms_on='msms', tol=(0.003, 0.005),
                       score_thd=0.0,
                       mode='w'):
        '''
        超大数据表内部的MS相似度计算
        分块推进
        之所以按块推进，是为了充分使用显卡的存储矩阵，若按行逐行向下计算，由于来回转储，时间开销大
        注意：返回的数据框也非常长，如果一致存储在内存，需要TB级的内存       
        性能：1200行数据，需要3.9469秒

        客户端在获得score长表后，最好按阈值过滤行，然后存储，以避免硬盘空间不够
        为方便后续的数据处理，客户端最好同时保存df和score，可以对应着过滤数据
        '''
        chunks = np.array_split(df.index, num_chunks)
        score_list = []
        n = len(chunks)
        f = open(csv_fname, mode=mode)
        f.write('row,col,score\n')
        f.flush()
        
        for i, chk in enumerate(chunks):
            # 计算当前块内部的相似度矩阵，转换为相似度长表
            que = df.loc[chk]
            match = cls.create_from_df(que, mz_on = mz_on, ms_on=ms_on)
            score = match.cosine_mx(tol=tol)
            score_df = match.to_long_df(que.index, que.index, score)
            score_df = score_df[(score_df['row'] < score_df['col']) &\
                                (score_df['score'] > score_thd)]
            score_df.to_csv(f, mode='a', header=False, index=False)
            f.flush()
            
            if i < n-1:
                #  计算当前块与其后所有块数据的相似度，转换为相似度长表
                # 下一个chunk的第一个索引值作为起点
                ref = df.loc[chunks[i+1][0]:]
                match = cls.create_from_df(ref, que,
                                               mz_on=mz_on,
                                               ms_on=ms_on)
                score = match.cosine_mx(tol=tol)
                score_df = match.to_long_df(ref.index, que.index, score)
                score_df = score_df[score_df['score'] > score_thd]
                score_df.to_csv(f, mode='a', header=False, index=False)
                f.flush()
            
        f.close()    
        


    @classmethod
    def create_from_df(cls, ref, que=None, 
                       mz_on     = 'precursormz',
                       ms_on     = 'msms',
                       ref_mz_on = None,
                       ref_ms_on = None,
                       que_mz_on = None,
                       que_ms_on = None,
                       sort_ms   = False
                      ):
        '''
        create objec from mzFrame
        '''
        if ref_mz_on is None:
            ref_mz_on = mz_on
        if ref_ms_on is None:
            ref_ms_on = ms_on
        if que_mz_on is None:
            que_mz_on = mz_on
        if que_ms_on is None:
            que_ms_on = ms_on
        if que is None:
            return cls(ref[ref_mz_on].values, ref[ref_ms_on].to_list(), sort_ms=sort_ms)
        else:
            return cls(ref[ref_mz_on].values, ref[ref_ms_on].to_list(),
                       que[que_mz_on].values, que[que_ms_on].to_list(),
                       sort_ms=sort_ms)

    def to_long(self, score_mx):
        '''
        transform score format
        '''
        rows, cols = np.indices(score_mx.shape)    
        return rows.ravel(), cols.ravel(), score_mx.ravel()

    def to_long_df(self, ref_index, que_index, score_mx):  
        """  
        通过 ref_index, query_index 和一个二维 scores 矩阵构建一个 DataFrame。  
        其中：  
            - ref_index 为行索引（一维，可视作形状 n）  
            - query_index 为列索引（一维，可视作形状 m）  
            - scores 为大小 n*m 的 2D 数组或矩阵  
    
        返回的 score_df 有三个列：  
            - ref_index: 对应原 ref_index 中的值  
            - query_index: 对应原 query_index 中的值  
            - score: scores 矩阵中的对应元素  
        """  
        rows, cols, scores = self.to_long(score_mx)

        score_df = pd.DataFrame({  
            "row": ref_index[rows].to_list(),  
            "col": que_index[cols].to_list(),  
            "score": np.round(scores, decimals=4)    
        })      
        return score_df 

    def to_csv(self, score_mx, filename, mode='w'):
        ''' 
        save similarity matrix to a csv file with three columns:
            row, col, value
        '''
        data = self.to_long(score_mx)        
        if mode.lower() in ("a", 'append'):  
            # 追加到已有文件末尾，不再写入列名  
            with open(filename, "ab") as f:  
                np.savetxt(  
                    f,  
                    data,  
                    delimiter=",",  
                    header="",      # 不再写入表头   
                    fmt=["%d", "%d", "%.4f"]  
                    )  
        else:  
            # 直接覆盖写入，并添加表头  
            with open(filename, "wb") as f:  
                np.savetxt(  
                    f,  
                    data,  
                    delimiter=",",  
                    header="row,col,value",  
                    fmt=["%d", "%d", "%.4f"]  
                )  

    @classmethod
    def warmup(cls):  
        """  
        在导入模块时自动调用的热身函数,用于预编译 _cosine_ 和 cosine_mx 函数  
        """  
        # 创建两个 GPUmatch 实例,并调用 cosine_mx 方法  
        ref_mz = np.array([100.0], dtype=np.float32)  
        ref_ms = [np.array([[100.0, 1000.0]], dtype=np.float32)]  
        que_mz = np.array([101.0], dtype=np.float32)  
        que_ms = [np.array([[101.0, 1000.0]], dtype=np.float32)]  

        gm1 = cls(ref_mz, ref_ms, que_mz, que_ms)  
        gm1.cosine_mx()  

        gm2 = cls(ref_mz, ref_ms)  
        gm2.cosine_mx()  

# 在导入模块时自动调用热身函数  
Match.warmup() 