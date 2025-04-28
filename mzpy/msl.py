# '''
# similarity中函数的并行加速版本
#     prepare_data 用于准备数据
#     函数名与ms包中保持一致
# '''
# from . import ms
# from . import similarity as sim
# from numba import njit, prange
# import numpy as np





# @njit
# def is_right_shape(que_list):
#     return que_list.ndim == 3 and que_list.shape[2] == 2 and \
#         que_list.shape[0] > 0 and que_list.shape[1] > 0  


# sim_get_scores = njit(sim.get_scores)

# @njit(parallel=True)  
# def get_scores(que_list, ref_list=None, tol=(0.003, 0.005), tol_type='Da'):
#     que_list = np.ascontiguousarray(que_list)
#     if ref_list is not None:
#         ref_list = np.ascontiguousarray(ref_list)

#     if not is_right_shape(que_list):
#         raise ValueError(f'Incompatible array shapes: {que_list.shape()}')
#     if ref_list is not None and not is_right_shape(ref_list):
#         raise ValueError(f'Incompatible array shapes: {ref_list.shape()}') 

#     val_len = 4  # 固定返回4个值  

#     x = que_list.shape[0]  

#     # 简单形状验证（Numba能支持）  
#     if que_list.ndim != 3 or que_list.shape[2] != 2:  
#         return np.empty((0, 2 + val_len))  

#     if ref_list is None:  
#         # 自身上三角  
#         row_count = x * (x + 1) // 2  
#         result = np.zeros((row_count, 2 + val_len), dtype=np.float64)  
#         idx = 0  
#         for i in prange(x):  
#             for j in range(i, x):  
#                 scores = sim_get_scores(que_list[i], que_list[j], tol)  
#                 # 这里假设scores为4元组  
#                 result[idx, 0] = i  
#                 result[idx, 1] = j  
#                 for k in range(val_len):  
#                     result[idx, 2 + k] = scores[k]  
#                 idx += 1  
#         return result  
#     else:  
#         # 计算所有对  
#         y = ref_list.shape[0]  
#         if ref_list.ndim != 3 or ref_list.shape[2] != 2:  
#             return np.empty((0, 2 + val_len))  
#         row_count = x * y  
#         result = np.zeros((row_count, 2 + val_len), dtype=np.float64)  
#         idx = 0  
#         for i in prange(x):  
#             for j in range(y):  
#                 scores = sim_get_scores(que_list[i], ref_list[j], tol)  
#                 for k in range(val_len):  
#                     result[idx, 2 + k] = scores[k]  
#                 result[idx, 0] = i  
#                 result[idx, 1] = j  
#                 idx += 1  
#         return result  







# @njit(cache=True)
# def _cosine_(que_sorted,
#              ref_sorted,
#              tol=(0.003, 0.005),
#              precursormz_compared=False,
#              penalty=True):
#     """计算余弦相似度的核心逻辑。
#     Assumption: que and ref's MS2 are already sorted in descending order by m/z, 
#                 and relative abundance has been applied.  
    
#     Args:  
#         que, MSList ndarray, shape (m, x, 2)  
#         ref, MSList ndarray, shape (n, y, 2)  
#         tol (tuple): 峰匹配的公差，包含两个值。  
#         to_sort,是否先排序再对齐mz。如果不对齐，计算的相似度不对
    
#     Returns:  
#         float: 计算得到的余弦相似度或 -1.0。  
#     """  
   
#     # 检查 precursor mz 差值
#     if precursormz_compared:
#         mz_diff = abs(que_sorted[0, 0] - ref_sorted[0, 0])  
#         if mz_diff > tol[0]:  
#             return 0.0

#     que_msms = que_sorted[1:]
#     ref_msms = ref_sorted[1:] 
#     que_msms = que_msms[que_msms[:, 0] != 0]
#     ref_msms = ref_msms[que_msms[:, 0] != 0]  

#     # MSList长度 
#     n_que = que_msms.shape[0]  
#     n_ref = ref_msms.shape[0]  

#     union_que = np.empty(n_que + n_ref, dtype=que_sorted.dtype)  
#     union_ref = np.empty(n_que + n_ref, dtype=ref_sorted.dtype)  

#     i = 0  # 从第二行开始  
#     j = 0  # 从第二行开始  
#     idx = 0  

#     while i < n_que and j < n_ref:  
#         mz_que, int_que = que_msms[i]  
#         mz_ref, int_ref = ref_msms[j]  

#         if abs(mz_que - mz_ref) <= tol[1]:  
#             union_que[idx] = int_que  
#             union_ref[idx] = int_ref  
#             i += 1  
#             j += 1  
#         elif mz_que > mz_ref:  
#             union_que[idx] = int_que  
#             union_ref[idx] = 0.0  
#             i += 1  
#         else:  
#             union_que[idx] = 0.0  
#             union_ref[idx] = int_ref  
#             j += 1  
#         idx += 1  

#     while i < n_que:  
#         union_que[idx] = que_msms[i, 1]  
#         union_ref[idx] = 0.0  
#         i += 1  
#         idx += 1  

#     while j < n_ref:  
#         union_que[idx] = 0.0  
#         union_ref[idx] = ref_msms[j, 1]  
#         j += 1  
#         idx += 1  

#     union_que = union_que[:idx]  
#     union_ref = union_ref[:idx] 
#     '''
#     union_que and union_ref are two arrays used to store the merged intensity values. 
#     Their specific meanings and roles are as follows:

#     - union_que: an array used to store the intensity values corresponding to 
#                     each matched mz in the query MS (que).
#         When a mz in `que` matches a mz in `ref`, `union_que` will store the intensity value in `que`.
    
#     The purpose of `union_que` and `union_ref` is to merge the intensity information of the two spectra 
#         into a unified structure for subsequent cosine similarity calculation.
#     Handling Unmatched Peaks: By setting the unmatched intensity values to 0, 
#         it ensures that the cosine similarity calculation is not affected.
#     '''
#     norm_que = np.sqrt(np.sum(union_que * union_que))  
#     norm_ref = np.sqrt(np.sum(union_ref * union_ref))  
#     if norm_que == 0 or norm_ref == 0:  
#         return 0.0  

#     dot_val = np.sum(union_que * union_ref)  
#     cos_val = dot_val / (norm_que * norm_ref)

#     if penalty:
#         lSpectrumCounter = np.sum(union_ref > 10)  

#         if lSpectrumCounter == 1:  
#             peakCountPenalty = 0.75  
#         elif lSpectrumCounter == 2:  
#             peakCountPenalty = 0.88  
#         elif lSpectrumCounter == 3:  
#             peakCountPenalty = 0.94  
#         elif lSpectrumCounter == 4:  
#             peakCountPenalty = 0.97  
#         else:  
#             peakCountPenalty = 1.0  

#         return cos_val * peakCountPenalty 
#     else:
#         return cos_val

# # @njit(cache=True, parallel=True)  
# # def _cosine_matrix_(que_msl,
# #                     ref_msl,
# #                     tol=(0.003, 0.005),
# #                     precursormz_compared=False,
# #                     penalty=True):  
# #     """  
# #     计算两个修整后的MSMS数组之间的MSMS相似度矩阵。  

# #     Args:  
# #         que_msl (ndarray): 查询的MSMS图谱数组，形状为(n, x, 2)。  
# #         ref_msl (ndarray): 参考的MSMS图谱数组，形状为(m, y, 2)。  
# #         tol (tuple): 峰匹配的公差，包含两个值。  
# #         precursormz_compared (bool): 是否比较前体质荷比。  

# #     Returns:  
# #         ndarray: MSMS相似度矩阵，形状为(n, m)。  
# #     """  
# #     num_que = que_msl.shape[0]  
# #     num_ref = ref_msl.shape[0]  
    
# #     # 初始化相似度矩阵  
# #     similarity_matrix = np.zeros((num_que, num_ref), dtype=np.float64)  

# #     for i in prange(num_que):  # 并行的触发点
# #         que = que_msl[i]
# #         for j in range(num_ref):                
# #             ref = ref_msl[j]  
# #             # 计算余弦相似度  
# #             similarity_matrix[i, j] = _cosine_(que, ref, tol, precursormz_compared, penalty)  

# #     return similarity_matrix 

# # @njit(cache=True, parallel=True)  
# # def _cosine_self_matrix_(que_msl, tol=(0.003, 0.005), precursormz_compared=False, penalty=True):  
# #     """  
# #     计算que_msl内部各个MSMS图谱之间的相似度矩阵。  

# #     Args:  
# #         que_msl (ndarray): 查询的MSMS图谱数组。  
# #         tol (tuple): 峰匹配的公差。  
# #         precursormz_compared (bool): 是否比较前体质荷比。  

# #     Returns:  
# #         ndarray: MSMS相似度矩阵，形状为(num_que, num_que)。  
# #     """  
# #     num_que = que_msl.shape[0]
# #     similarity_matrix = np.zeros((num_que, num_que), dtype=np.float64)  

# #     for i in prange(num_que):  
# #         similarity_matrix[i, i] = 1.0  # 自身相似度为1.0  
# #         que_1 = que_msl[i]  
# #         for j in range(i + 1, num_que):              
# #             que_2 = que_msl[j]  
# #             similarity = _cosine_(que_1, que_2, tol, precursormz_compared, penalty)  
# #             similarity_matrix[i, j] = similarity  
# #             similarity_matrix[j, i] = similarity  # 复制到对称位置  

# #     return similarity_matrix


# # class MSList(np.ndarray):
# #     def __new__(cls, ms_arr):
# #         '''
# #         Organize the nested MSMS spectrum array into a NumPy 3D array with shape (n, x, 2).
# #             The x dimension is determined by the number of rows in the MSMS spectrum with the most fragment ion.
# #             And other MS are padded with zeros to complete the data.
# #         param:
# #             msms_arr: a list or NumPy array, where each element is a 2D array of MSMS.
# #                         The length of each MSMS may vary.
# #         '''
# #         # 填充补0，MSMS数组等长度
# #         max_length = 0
# #         for i, ms in enumerate(ms_arr):
# #             if not isinstance(ms, ms.MSdata):
# #                 ms_arr[i] = ms.asMSdata(ms)
# #             if ms.shape[0] > max_length:
# #                 max_length = ms.shape[0]


        
# #         # 创建新的 NumPy 数组，填充零
# #         n = len(ms_arr)
# #         adjusted_arr = np.zeros((n, max_length, 2), dtype=np.float64)
        
# #         # 填充数组
# #         for i in range(n):
# #             precursor = msms_arr[i][0]
# #             msms = msms_arr[i][1:]
# #             msms = sorted(msms, key=lambda x: x[0], reverse=True)
# #             # msms_sorted =  msms[np.argsort(msms[:, 0])[::-1]]
# #             if to_RA:
# #                 msms = np.asarray(msms)
# #                 msms = to_RA(msms) 

# #             data = np.vstack((precursor, msms))
# #             current_length = len(msms_arr[i])
# #             adjusted_arr[i, :current_length] = np.asarray(data)       

# #         return np.ascontiguousarray(adjusted_arr).view(cls)      
    
# #     def __array_finalize__(self, obj):
# #         # 这个方法处理切片等操作
# #         if obj is None: return
    
# #     @classmethod
# #     def create(cls, mz_list, msms_list, intensity=100):
# #         '''
# #         创建 MSList2 对象，每个 MSMS 谱图的头部添加母离子质荷比和丰度（默认为100）
        
# #         参数:
# #             precursor_mzs: 浮点数数组，记录母离子的质荷比
# #             msms_list: MSMS图谱的列表，每个元素是一个二维数组
# #         '''
# #         # 确保两个输入参数长度相等
# #         if (not isinstance(mz_list, (list, np.ndarray))) or \
# #            (not isinstance(msms_list, (list, np.ndarray))):
# #             raise TypeError(f'mz_list must be a list or numpy array.')
# #         if len(mz_list) != len(msms_list):
# #             raise ValueError("precursor_mzs 和 msms_list 的长度必须相等")
        
# #         ms = []
# #         for i, mz in enumerate(mz_list):
# #             ms.append(np.vstack(([mz, intensity], msms_list[i])))
        
# #         return cls(ms)
    
# #     def to_RA(self):
# #         '''
# #         transfer intensity data to ralative abundance
# #         RA: relative abundance
# #         '''
# #         # data = []
# #         msms = to_RA(self[1:])
# #         return np.vstack((self[0], msms))

                
# #     def sort_ms2(self, ascending=False):
# #         '''
# #         sort each MSMS according mz
# #         '''
# #         pass

# #     def match(self,
# #               ref_msl=None,
# #               tol=(0.003, 0.005),
# #               precursormz_compared=False,
# #               penalty=True):
# #         if ref_msl is None:
# #             return _cosine_self_matrix_(self,
# #                                         tol=tol,
# #                                         precursormz_compared=precursormz_compared,
# #                                         penalty=penalty)
# #         else:
# #             if not isinstance(ref_msl, self.__class__):
# #                 raise TypeError(f'ref_msl nust be an object of {self.__class__.__name__}')
# #             return _cosine_matrix_(self,
# #                                 ref_msl,
# #                                 tol=tol,
# #                                 precursormz_compared=precursormz_compared,
# #                                 penalty=penalty)
    
# #     def warmup(self, size=10):
# #         import time
# #         n = min(size, self.shape[0])
# #         if n > 0:
# #             # 从中随机抽取2个成员  
# #             random_indices = np.random.choice(self.shape[0], size=2, replace=False)  
# #             random_samples = self[random_indices] 
# #             start_time = time.time()
# #             a = _cosine_matrix_(random_samples,
# #                             random_samples)
            
# #             b = _cosine_self_matrix_(random_samples)
# #             end_time = time.time() 
# #             print(f'time used: {end_time - start_time} sec.')           


      

# # ############    preheat    ############

# # def warmup_numba():
# #     # 创建示例数据
# #     example_que_msl = np.random.rand(2, 3, 2)  # 示例查询 MSList
# #     example_ref_msl = np.random.rand(3, 3, 2)  # 示例参考 MSList
# #     example_tol = (0.003, 0.005)
# #     example_precursormz_compared = True
# #     example_penalty = True
# #     print('warm up numba ...', end=' ')
# #     # 调用 _cosine_matrix_ 和 _cosine_self_matrix_ 来触发编译
# #     import time
# #     start_time = time.time()
# #     _cosine_matrix_(example_que_msl,
# #                     example_ref_msl,
# #                     example_tol,
# #                     example_precursormz_compared,
# #                     example_penalty)
    
# #     _cosine_self_matrix_(example_que_msl,
# #                          example_tol,
# #                          example_precursormz_compared,
# #                          example_penalty)
# #     end_time = time.time()
# #     print(f'time used: {end_time - start_time} sec.', end='\t')
# #     print('successfully.')

# # # 在程序启动时调用预热函数
# # warmup_numba()
# ###################################