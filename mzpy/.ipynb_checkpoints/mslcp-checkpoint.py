import cupy as cp
import numpy as np
from .msl import MSList


def cosine_batch(que_batch, ref_batch, tol=(0.003, 0.005), precursormz_compared=True):  
    """批量计算余弦相似度  
    
    Args:  
        que_batch (list of cupy.ndarray): 查询质谱数据批次  
        ref_batch (list of cupy.ndarray): 参考质谱数据批次  
        tol (tuple): 质荷比匹配公差  
        
    Returns:  
        cupy.ndarray: 相似度矩阵，形状为(len(que_batch), len(ref_batch))  
    """  
    batch_que_size = len(que_batch)  
    batch_ref_size = len(ref_batch)  
    
    # 预分配结果矩阵  
    result_matrix = cp.full((batch_que_size, batch_ref_size), -1.0)  
    
    # 批量处理precursor m/z筛选  
    if precursormz_compared:  
        que_precursor_mzs = cp.array([spec[0, 0] for spec in que_batch])  
        ref_precursor_mzs = cp.array([spec[0, 0] for spec in ref_batch])  
        
        # 计算所有组合的m/z差异  
        mz_diff_matrix = cp.abs(que_precursor_mzs[:, cp.newaxis] - ref_precursor_mzs[cp.newaxis, :])  
        valid_pairs = mz_diff_matrix <= tol[0]  
    else:  
        valid_pairs = cp.ones((batch_que_size, batch_ref_size), dtype=bool)  
    
    # 对有效对进行余弦相似度计算  
    for i in range(batch_que_size):  
        for j in range(batch_ref_size):  
            if valid_pairs[i, j]:  
                result_matrix[i, j] = _cosine_single_pair(que_batch[i], ref_batch[j], tol)  
    
    return result_matrix  


def _cosine_single_pair(que_sorted, ref_sorted, tol=(0.003, 0.005)):  
    """计算余弦相似度的核心逻辑。  
    Assumption: que_sorted and ref_sorted's MS2 are already sorted in descending order by m/z,   
                and relative abundance has been applied.  
    
    Args:  
        que_sorted (cupy.ndarray): mz从大到小排序后的查询数据。  
        ref_sorted (cupy.ndarray): mz从大到小排序后的参考数据。  
        tol (tuple): 峰匹配的公差，包含两个值。  
    
    Returns:  
        float matrix: similarity matrix
    """  
    # 计算 MS2 的余弦相似度  
    n_que = que_sorted.shape[0]  
    n_ref = ref_sorted.shape[0]  

    union_que = cp.empty(n_que + n_ref, dtype=que_sorted.dtype)  
    union_ref = cp.empty(n_que + n_ref, dtype=ref_sorted.dtype)  

    i = 1  # 从第二行开始  
    j = 1  # 从第二行开始  
    idx = 0  

    while i < n_que and j < n_ref:  
        mz_que, int_que = que_sorted[i]  
        mz_ref, int_ref = ref_sorted[j]  
        if cp.abs(mz_que - mz_ref) <= tol[1]:  
            union_que[idx] = int_que  
            union_ref[idx] = int_ref  
            i += 1  
            j += 1  
        elif mz_que > mz_ref:  
            union_que[idx] = int_que  
            union_ref[idx] = 0.0  
            i += 1  
        else:  
            union_que[idx] = 0.0  
            union_ref[idx] = int_ref  
            j += 1  
        idx += 1  

    while i < n_que:  
        union_que[idx] = que_sorted[i, 1]  
        union_ref[idx] = 0.0  
        i += 1  
        idx += 1  

    while j < n_ref:  
        union_que[idx] = 0.0  
        union_ref[idx] = ref_sorted[j, 1]  
        j += 1  
        idx += 1  

    union_que = union_que[:idx]  
    union_ref = union_ref[:idx]   

    norm_que = cp.sqrt(cp.sum(union_que * union_que))  
    norm_ref = cp.sqrt(cp.sum(union_ref * union_ref))  
    if norm_que == 0 or norm_ref == 0:  
        return 0.0  

    dot_val = cp.sum(union_que * union_ref)  
    cos_val = dot_val / (norm_que * norm_ref)  

    lSpectrumCounter = cp.sum(union_ref > 10)  

    if lSpectrumCounter == 1:  
        peakCountPenalty = 0.75  
    elif lSpectrumCounter == 2:  
        peakCountPenalty = 0.88  
    elif lSpectrumCounter == 3:  
        peakCountPenalty = 0.94  
    elif lSpectrumCounter == 4:  
        peakCountPenalty = 0.97  
    else:  
        peakCountPenalty = 1.0  

    return cos_val * peakCountPenalty



# class MSList_cp():  
#     def __init__(self, precursormz_list, msms_list, precursor_intensity=100):  
#         """  
#         初始化 CuPyMSList 实例，创建一个包含母离子和对应子离子信息的数组。  
        
#         参数:  
#         precursormz_list: 母离子的质荷比列表 (list of float)  
#         msms_list: 对应的子离子信息列表，每个元素为二维数组 (list of cupy.ndarray)  
#         precursor_intensity: 母离子的强度，默认为 100  
#         """  
#         msl = MSList(precursormz_list,
#                      msms_list,
#                      precursor_intensity=precursor_intensity) 
#         self.data = [cp.asarray(arr) for arr in msl.data] 
        
#     def compute_similarity(self, ref_ms_list, tol=(0.003, 0.005), precursormz_compared=True):  
#         '''  
#         计算 self 列表与 ref_ms_list 之间的余弦相似度  
#         ref_ms_list, an instance of CuPyMSList  
#         return:  
#             a matrix of cosine similarity  
#         '''  
#         if not isinstance(ref_ms_list, self.__class__):  
#             raise TypeError(f'ref_ms_list must be an instance of {self.__class__.__name__}')  
        
#         n_que = self.data.shape[0]  
#         n_ref = ref_ms_list.data.shape[0]  
#         similarity_matrix = cp.full((n_que, n_ref), -1.0)  # 初始化相似度矩阵，默认值为 -1.0  

#         for i in range(n_que):  
#             for j in range(n_ref):  
#                 similarity_matrix[i, j] = _cosine_(self.data[i],  
#                                                    ref_ms_list.data[j],  
#                                                    tol,
#                                                    precursormz_compared)  

#         return similarity_matrix  

#     def compute_similarity_self(self, tol=(0.003, 0.005), precursormz_compared=True):  
#         '''  
#         自相似度计算  
#         '''       
#         n = len(self.data)
#         similarity_matrix = cp.full((n, n), -1.0)  # 初始化相似度矩阵，默认值为 -1.0  
#         for i in range(n):  
#             for j in range(n):   
#                 if i < j:  # 上三角  
#                     similarity_matrix[i, j] = _cosine_(self.data[i],  
#                                                        self.data[j],  
#                                                        tol,
#                                                        precursormz_compared)   
#                 elif i == j:  # 对角线  
#                     similarity_matrix[i, j] = 1.0  
#                 elif i > j:  
#                     similarity_matrix[i, j] = similarity_matrix[j, i]  
#         return similarity_matrix  

#     def to_cupy(self):  
#         """将数据转换为 CuPy 数组"""  
#         return cp.ascontiguousarray(self.data)    


################    preheat    ################         

# # 预热代码  
# def preheat_cupy(): 
#     print('preheating _cosine_ function')  
#     # 创建一个简单的CuPy数组并执行一些操作  
#     precursormz_list = [100.0, 150.0]  # 示例质荷比列表  
#     msms_list = [[[100.0, 50.0], [101.0, 25.0]],  # 示例子离子信息  
#                  [[150.0, 75.0], [151.0, 30.0]]]  
    
#     # 创建MSList_cp的实例  
#     ms_list = MSList_cp(precursormz_list, msms_list)  

# # 调用预热函数  
# preheat_cupy()  
