'''
Function Family for Calculating Mass Spectral Similarity

ref:
    MsdialWorkbench/src/Common/CommonStandard/Algorithm/Scoring/MsScanMatching.cs
    https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/MSDIAL5/MsdialCore/Algorithm/Annotation/Ms2MatchCalculator.cs

'''
from numba import njit, prange
import numpy as np

from .ms import MSdata

__WARMED__ = False # 用于标识是否执行过预热函数



#### data prepare functions
 
def join(precursormz, msms, intensity=100, to_normlized=True):  
    """  
    拼接母离子和碎片离子矩阵形成一组MS  
    """  
    data = np.vstack(([[precursormz, intensity]], msms))  
    return MSdata(data, to_normalized=to_normlized)  
 
def join_array(mz_list, msms_list, intensity=100, to_normlized=True):  
    """  
    拼接mz数组和msms数组形成MS数组  
    """  
    if len(mz_list) != len(msms_list):  
        raise ValueError(f'Inconsistent Length: {len(mz_list)} vs. {len(msms_list)}')  

    data = [join(mz, msms_list[i], intensity, to_normlized)   
            for i, mz in enumerate(mz_list)]  
    return data  


def pad_array(ms_arr):  
    """  
    用0补齐列表中的所有MSdata数据数据，并形成三维数组  
    """  
    msdata_arrays = [  
        MSdata(item) if isinstance(item, (list, np.ndarray)) else item  
        for item in ms_arr  
    ]  
    lengths = [len(msdata) for msdata in msdata_arrays]  
    max_length = max(lengths)  

    padded_array = np.empty((len(msdata_arrays), max_length, 2))  

    for i, msdata in enumerate(msdata_arrays):  
        padded_array[i, :len(msdata), :] = msdata  

    return padded_array  


def prepare_ms_list(mz_list, msms_list, intensity=100, to_normlized=True):  
    data = join_array(mz_list, msms_list, intensity, to_normlized)  
    return pad_array(data)  

#######################################################################################

## similarity caculator        

def warmup(que_list=None, ref_list=None):  
    """pre heat numba functions""" 
    # 模拟简单小谱用于预热 
    global __WARMED__
       
    if not __WARMED__:
        print('preheating numba functions ... ', end='\t') 
        if que_list is None:
            que = np.array([[100.0, 200.0], [50.0, 80.0], [60.0, 90.0]], dtype=np.float64)  
            que_list = np.stack([que, que * 1.01], axis=0)
        if ref_list is None:
            ref = np.array([[100.0, 210.0], [50.0, 85.0], [61.0, 95.0]], dtype=np.float64)  
            ref_list = np.stack([ref, ref * 0.99], axis=0)  
        tol = (0.003, 0.005)  

        # 调用search_matched_peaks和get_scores  
        # _ = get_union_peaks(que, ref, tol)  
        # _ = get_scores(que, ref, tol)      

        _ = get_scores_batch(que_list, ref_list, tol)  
        # _ = get_scores_batch(que_list, None, tol)  
        print('accomplished.')
        __WARMED__ = True


@njit
def get_union_peaks(que_sorted, ref_sorted, tol=(0.003, 0.005)):  
    """  
    合并并匹配两条已按 m/z 降序排序的 MS2 光谱，返回 shape=(N, 2) 的强度矩阵。  
    """  
    if abs(que_sorted[0, 0] - ref_sorted[0, 0]) > tol[0]:  
        return np.empty((0, 2), dtype=que_sorted.dtype)  

    que = que_sorted[que_sorted[:, 0] > np.float32(0.0)]  
    ref = ref_sorted[ref_sorted[:, 0] > np.float32(0.0)]  

    n_que = que.shape[0]  
    n_ref = ref.shape[0]  
    max_len = n_que + n_ref - 2  # 不计母离子，最大可能匹配数  

    union = np.zeros((max_len, 2), dtype=que.dtype)  
    i, j, idx = 1, 1, 0  ## i, j计数器从1开始，跳过母离子

    while i < n_que and j < n_ref:  
        mz_q, int_q = que[i]  
        mz_r, int_r = ref[j]  
        if abs(mz_q - mz_r) <= tol[1]:  
            union[idx, 0] = int_q  
            union[idx, 1] = int_r  
            i += 1  
            j += 1  
        elif mz_q > mz_r:  
            union[idx, 0] = int_q  
            union[idx, 1] = 0.0  
            i += 1  
        else:  
            union[idx, 0] = 0.0  
            union[idx, 1] = int_r  
            j += 1  
        idx += 1  

    while i < n_que:  
        union[idx, 0] = que[i, 1]  
        union[idx, 1] = 0.0  
        i += 1  
        idx += 1  

    while j < n_ref:  
        union[idx, 0] = 0.0  
        union[idx, 1] = ref[j, 1]  
        j += 1  
        idx += 1  

    return union[:idx]  


@njit
def get_scores(que, ref, tol=(0.003, 0.005)):
    '''
    文大模型：这段代码如何避免除以0错误
    '''

    arr = get_union_peaks(que, ref, tol)
    matched_mask = (arr[:, 0] != np.float32(0)) & (arr[:, 1] != np.float32(0))
    matched_count = np.sum(matched_mask)

    if matched_count == 0:
        return 0, 0.0, 0.0
    # 匹配峰的点积

    matched_product = np.where(  
                        matched_mask, # 都非零  
                        arr[:, 0] * arr[:, 1], 0
                    )  
    matched_product_total = matched_product.sum()

    # 不匹配峰的点积
    unmatched_product = np.where(  
                    (arr[:, 0] == 0) | (arr[:, 1] == 0), # 都非零  
                    arr[:, 0]**2 + arr[:, 1]**2, 0
                ) 
    bonanza_denominator = matched_product_total + unmatched_product.sum()  
    if bonanza_denominator == 0:  
        bonanza = 0.0  
    else:  
        bonanza = matched_product_total / bonanza_denominator  
    
    mask = (arr[:, 0] != 0) & (arr[:, 1] != 0)  
    result1 = np.sqrt(np.sum(arr[mask, 0] ** 2)) 
    result2 = np.sqrt(np.sum(arr[mask, 1] ** 2)) 
    
    # 归一化的余弦相似度
    denominator_cosine = result1 * result2  
    if denominator_cosine == 0:  
        cosine = 0.0  
    else:  
        cosine = matched_product_total / denominator_cosine  

    return matched_count, bonanza, cosine



@njit(parallel=True)  
def get_scores_batch(que_list, ref_list=None, tol=(0.003, 0.005)):  
    n1 = que_list.shape[0]  

    # ref_list为空时，行列都对应que_list  
    if (ref_list is None) or (ref_list.shape[0] == 0):  
        matched_mx = np.zeros((n1, n1), dtype=np.int32)  
        bonanza_mx = np.zeros((n1, n1), dtype=np.float32)  
        cosine_mx  = np.zeros((n1, n1), dtype=np.float32)  

        for i in prange(n1):  
            for j in range(i, n1):  
                if i == j:
                    matched_count =  np.sum(que_list[i][:, 0] > 0)  
                    bonanza = 1.0  
                    cosine = 1.0  
                else:  
                    matched_count, bonanza, cosine = get_scores(que_list[i],
                                                                que_list[j],
                                                                tol) 
                     
                matched_mx[i, j] = matched_count  
                bonanza_mx[i, j] = bonanza  
                cosine_mx[i, j] = cosine 

    else:  
        n2 = ref_list.shape[0]  
        matched_mx = np.zeros((n1, n2), dtype=np.int32)  
        bonanza_mx = np.zeros((n1, n2), dtype=np.float32)  
        cosine_mx  = np.zeros((n1, n2), dtype=np.float32)  

        for i in prange(n1):  
            for j in range(n2):  
                matched_count, bonanza, cosine = get_scores(que_list[i],
                                                            ref_list[j],
                                                            tol)  
                matched_mx[i, j] = matched_count  
                bonanza_mx[i, j] = bonanza  
                cosine_mx[i, j] = cosine  
    return matched_mx, bonanza_mx, cosine_mx


 