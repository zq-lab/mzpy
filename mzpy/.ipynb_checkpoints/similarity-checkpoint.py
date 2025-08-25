'''
Function Family for Calculating Mass Spectral Similarity

ref:
    MsdialWorkbench/src/Common/CommonStandard/Algorithm/Scoring/MsScanMatching.cs
    https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/MSDIAL5/MsdialCore/Algorithm/Annotation/Ms2MatchCalculator.cs

'''

import numba
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

# def warmup(que_list=None, ref_list=None):  
#     """pre heat numba functions""" 
#     # 模拟简单小谱用于预热 
#     global __WARMED__
       
#     if not __WARMED__:
#         print('preheating numba functions ... ', end='\t') 
#         if que_list is None:
#             que = np.array([[100.0, 200.0], [50.0, 80.0], [60.0, 90.0]], dtype=np.float64)  
#             que_list = np.stack([que, que * 1.01], axis=0)
#         if ref_list is None:
#             ref = np.array([[100.0, 210.0], [50.0, 85.0], [61.0, 95.0]], dtype=np.float64)  
#             ref_list = np.stack([ref, ref * 0.99], axis=0)  
#         tol = (0.003, 0.005)  

#         # 调用search_matched_peaks和get_scores  
#         # _ = get_union_peaks(que, ref, tol)  
#         # _ = get_scores(que, ref, tol)      

#         _ = get_scores_batch(que_list, ref_list, tol)  
#         # _ = get_scores_batch(que_list, None, tol)  
#         print('accomplished.')
#         __WARMED__ = True


@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_bonanza_score(union_peaks):
    """
    计算两个质谱之间的Bonanza得分
    
    参数:
        union_peaks: 由get_union_peaks函数返回的形状为(n, 2)的numpy数组
                    表示两个谱图的匹配峰和不匹配峰信息

    计算公式：
        匹配峰的强度乘积和 / (匹配峰强度乘积和 + 未匹配峰的强度平方和1 + 未匹配峰的强度平方和2)
                    
    返回值:
        Bonanza得分: 范围[0,1]，其中1表示完全相似，0表示完全不相似
    """
    if union_peaks.shape[0] == 0:
        return 0.0
    
    left_intensity = union_peaks[:, 0]
    right_intensity = union_peaks[:, 1]
    
    # 找出匹配峰和不匹配峰的索引
    matched_indices = (left_intensity > 0) & (right_intensity > 0)
    left_only_indices = (left_intensity > 0) & (right_intensity == 0)
    right_only_indices = (left_intensity == 0) & (right_intensity > 0)
    
    # 如果没有匹配峰，得分为0
    if np.sum(matched_indices) == 0:
        return 0.0
    
    # 计算匹配峰强度乘积和
    matched_product_sum = np.sum(left_intensity[matched_indices] * right_intensity[matched_indices])
    
    # 计算左谱图中未匹配峰的强度平方和
    left_unmatched_sum = np.sum(left_intensity[left_only_indices] ** 2)
    
    # 计算右谱图中未匹配峰的强度平方和
    right_unmatched_sum = np.sum(right_intensity[right_only_indices] ** 2)
    
    # 计算Bonanza得分
    denominator = matched_product_sum + left_unmatched_sum + right_unmatched_sum
    
    # 避免除零错误
    if denominator == 0:
        return 0.0
    
    score = matched_product_sum / denominator
    
    # 确保结果在[0, 1]范围内
    return max(0.0, min(1.0, score))


@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_matched_num(union_peaks):
    '''
    返回匹配峰的个数
    '''
    return np.sum((union_peaks[:, 0] > 0) & (union_peaks[:, 1] > 0))


@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_matched_peaks_score(union_peaks):
    """
    计算两个质谱之间的峰匹配得分
    
    参数:
        union_peaks: 由get_union_peaks函数返回的形状为(n, 2)的numpy数组
                     表示两个谱图的匹配峰和不匹配峰信息

    计算过程:
        计算ms_left中有多少峰能在ms_right中找到匹配，返回匹配百分比
        结果为[0,1]之间的值，0表示没有相似性，1表示完全一致
        返回两个值：[0]匹配峰的比例，[1]匹配峰的数量
                     
    返回值:
        匹配峰比例: 匹配峰数量占总峰数的比例 [0,1]
    """
    if union_peaks.shape[0] == 0:
        return 0.0
    
    # 计算匹配峰数量 (两列都大于0的行数)
    matched_peaks = np.sum((union_peaks[:, 0] > 0) & (union_peaks[:, 1] > 0))
    
    # 计算左谱图中的峰数量
    left_peaks = np.sum(union_peaks[:, 0] > 0)
    
    # 计算右谱图中的峰数量
    right_peaks = np.sum(union_peaks[:, 1] > 0)
    
    # 根据C#代码逻辑，取两者中的最大值作为分母
    denominator = max(left_peaks, right_peaks)
    
    # 避免除零错误
    if denominator == 0:
        return 0.0
    
    # 计算匹配比例
    match_ratio = float(matched_peaks) / denominator
    
    return match_ratio

@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_modified_dot_product_score(union_peaks):
    """
    计算两个质谱之间的修正点积得分
    
    参数:
        union_peaks: 由get_union_peaks函数返回的形状为(n, 2)的numpy数组
                    表示两个谱图的匹配峰和不匹配峰信息
    
    计算过程:
        查找两个谱图之间的匹配峰
        计算匹配峰的强度乘积之和
        计算两个谱图峰值强度的平方和作为标准化因子
        最终得分 = 乘积和 / (标准化因子1的平方根 * 标准化因子2的平方根)
                    
    返回值:
        修正点积得分: 范围[0,1]，其中1表示完全相似，0表示完全不相似
    """
    if union_peaks.shape[0] == 0:
        return 0.0
    
    left_intensity = union_peaks[:, 0]
    right_intensity = union_peaks[:, 1]
    
    # 找出匹配峰的索引（两个谱图中都有值的峰）
    matched_indices = (left_intensity > 0) & (right_intensity > 0)
    
    # 如果没有匹配峰，得分为0
    if np.sum(matched_indices) == 0:
        return 0.0
    
    # 计算匹配峰强度乘积和
    dot_product = np.sum(left_intensity[matched_indices] * right_intensity[matched_indices])
    
    # 计算标准化因子
    norm1 = np.sqrt(np.sum(left_intensity ** 2))
    norm2 = np.sqrt(np.sum(right_intensity ** 2))
    
    # 避免除零错误
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    # 计算修正点积得分
    score = dot_product / (norm1 * norm2)
    
    # 确保结果在[0, 1]范围内
    return max(0.0, min(1.0, score))

@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_simple_dot_product(union_peaks):
    """
    计算两个质谱之间的简单点积相似度,也被称为余弦相似度
    
    参数:
        union_peaks: 由get_union_peaks函数返回的形状为(n, 2)的numpy数组
                    表示两个谱图的匹配峰和不匹配峰信息
                    
    返回值:
        简单点积相似度: 范围[0,1]，其中1表示完全相似，0表示完全不相似
    """
    if union_peaks.shape[0] == 0:
        return 0.0
    
    # 提取强度值
    left_intensity = union_peaks[:, 0]
    right_intensity = union_peaks[:, 1]
    
    # 计算点积
    dot_product = np.sum(left_intensity * right_intensity)
    
    # 计算标准化因子
    norm_left = np.sqrt(np.sum(left_intensity ** 2))
    norm_right = np.sqrt(np.sum(right_intensity ** 2))
    
    # 避免除零错误
    if norm_left == 0 or norm_right == 0:
        return 0.0
    
    # 计算标准化点积（余弦相似度）
    score = dot_product / (norm_left * norm_right)
    
    # 确保结果在[0, 1]范围内
    return max(0.0, min(1.0, score))

@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_reverse_dot_product(union_peaks):
    """
    计算逻辑不对，相同的质谱的相似度也只有0.067
    计算两个质谱之间的反向点积相似度
    用于比较实验谱与标准谱
    基于Stein等人的方法（1999年发表）
    
    参数:
        union_peaks: 由get_union_peaks函数返回的形状为(n, 2)的numpy数组
                   表示两个谱图的匹配峰和不匹配峰信息

    计算过程:    
        包括质量加权和强度归一化
        计算公式: (协方差的平方) / (标量M * 标量R * 峰数惩罚因子)
                   
    返回值:
        反向点积相似度: 范围[0,1]，其中1表示完全相似，0表示完全不相似
    """
    if union_peaks.shape[0] == 0:
        return 0.0
    
    # 提取强度值
    left_intensity = union_peaks[:, 0]
    right_intensity = union_peaks[:, 1]
    
    # 找出匹配峰的索引
    matched_indices = (left_intensity > 0) & (right_intensity > 0)
    
    if np.sum(matched_indices) == 0:
        return 0.0
    
    # 1. 归一化强度值
    sum_left = np.sum(left_intensity)
    sum_right = np.sum(right_intensity)
    
    if sum_left == 0 or sum_right == 0:
        return 0.0
    
    norm_left = left_intensity / sum_left
    norm_right = right_intensity / sum_right
    
    # 2. 计算协方差（点积）
    covariance = np.sum(norm_left[matched_indices] * norm_right[matched_indices])
    
    # 3. 计算归一化因子
    scalar_m = np.sqrt(np.sum(norm_left ** 2))
    scalar_r = np.sqrt(np.sum(norm_right ** 2))
    
    if scalar_m == 0 or scalar_r == 0:
        return 0.0
    
    # 4. 使用get_matched_peaks_score函数计算峰数惩罚因子
    # 获取匹配峰比例作为峰数惩罚因子
    # peak_penalty = get_matched_peaks_score(union_peaks)
    
    # 5. 计算最终得分
    # 按照公式: (covariance * covariance) / (scalarM * scalarR * peakPenaltyFactor)
    # score = (covariance * covariance) / (scalar_m * scalar_r * peak_penalty) if peak_penalty > 0 else 0.0
    score = (covariance * covariance) / (scalar_m * scalar_r)
    # 确保结果在[0, 1]范围内
    return max(0.0, min(1.0, score))

@njit(numba.float64(numba.float64[:,:]), cache=True, fastmath=True, nogil=True)
def get_entropy_similarity(union_peaks):
    """
    计算两个质谱之间的谱熵相似度
    
    参数:
        union_peaks: 由get_union_peaks函数返回的形状为(n, 2)的numpy数组
                    表示两个谱图的匹配峰和不匹配峰信息

    计算过程:
        计算两个谱图的组合谱图熵
        分别计算每个谱图的熵
        通过公式: 1 - (2 * entropy12 - entropy1 - entropy2) * 0.5 计算相似度
                    
    返回值:
        谱熵相似度: 范围[0,1]，其中1表示完全相似，0表示完全不相似
    """
    if union_peaks.shape[0] == 0:
        return 0.0
    
    left_intensity = union_peaks[:, 0]
    right_intensity = union_peaks[:, 1]
    
    # 检查两个谱图是否有有效强度
    sum_left = np.sum(left_intensity)
    sum_right = np.sum(right_intensity)
    if sum_left == 0 or sum_right == 0:
        return 0.0
    
    # 计算组合谱图的强度
    combined_intensity = left_intensity + right_intensity
    
    # 归一化为概率分布
    p_left = left_intensity / sum_left
    p_right = right_intensity / sum_right
    p_combined = combined_intensity / np.sum(combined_intensity)
    
    # 计算熵，避免log(0)的情况
    def _entropy(p):
        # 只计算非零概率的熵
        p_non_zero = p[p > 0]
        if len(p_non_zero) == 0:
            return 0.0
        return -np.sum(p_non_zero * np.log(p_non_zero))
    
    entropy1 = _entropy(p_left)
    entropy2 = _entropy(p_right)
    entropy12 = _entropy(p_combined)
    
    # 计算谱熵相似度
    similarity = 1.0 - 0.5 * (2 * entropy12 - entropy1 - entropy2)
    
    # 确保结果在[0, 1]范围内
    return max(0.0, min(1.0, similarity))


def get_scores(ms_left,
               ms_right,
               tol=(0.003, 0.005),
               weights=None,
               method="weighted_average", 
               include_scores=None,
               exclude_scores=None):
    """
    基于多种相似度指标计算多种相似度
    
    参数:
        ms_left: 左谱图数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
        ms_right: 右谱图数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
        tol: 峰匹配的质荷比容差，默认为0.005
        weights: 可选，字典类型，指定每种相似度指标的权重
                如果为None，则使用默认权重
        method: 可选，字符串类型，指定综合打分的方法
                "weighted_average": 加权平均
                "geometric_mean": 几何平均(所有得分的乘积开n次方)
                "harmonic_mean": 调和平均
                "max": 取最大值
                "min": 取最小值
        include_scores: 可选，列表类型，仅包含指定的相似度指标
        exclude_scores: 可选，列表类型，排除指定的相似度指标
                
    返回值:
        scores_dict: 字典类型，包含各个单独的相似度得分和综合得分("total_score")
    """
    # 获取合并峰和质荷比信息
    union_peaks = get_union_ms(ms_left, ms_right, tol[0], tol[1])
    
    # 初始化可用的相似度计算函数
    available_scores = {
        "matched_peaks_ratio": lambda: get_matched_peaks_score(union_peaks),
        "spectral_entropy": lambda: get_entropy_similarity(union_peaks),
        "modified_dot_product": lambda: get_modified_dot_product_score(union_peaks),
        "bonanza_score": lambda: get_bonanza_score(union_peaks),
        "reverse_dot_product": lambda: get_reverse_dot_product(union_peaks),
        "simple_dot_product": lambda: get_simple_dot_product(union_peaks)
    }
    
    # 设置默认权重(如果未提供权重)
    if weights is None:
        weights = {
            "matched_peaks_ratio": 1.0,
            "spectral_entropy": 1.0,
            "modified_dot_product": 1.5,
            "bonanza_score": 1.5,
            "reverse_dot_product": 1.0,
            "simple_dot_product": 0.75
        }
    
    # 处理包含和排除列表
    score_keys = set(available_scores.keys())
    
    if include_scores is not None:
        score_keys = score_keys.intersection(set(include_scores))
    
    if exclude_scores is not None:
        score_keys = score_keys.difference(set(exclude_scores))
    
    # 计算指定的相似度得分
    scores = {}
    for score_name in score_keys:
        if score_name in available_scores:
            scores[score_name] = available_scores[score_name]()
    
    # 如果没有计算任何得分，返回空字典
    if not scores:
        scores["total_score"] = 0.0
        return scores
    
    # 根据指定方法计算综合得分
    if method == "weighted_average":
        # 计算总权重
        total_weight = sum(weights.get(name, 1.0) for name in scores)
        
        # 避免除零错误
        if total_weight == 0:
            scores["total_score"] = 0.0
            return scores
        
        # 计算加权平均
        composite_score = sum(score * weights.get(name, 1.0) for name, score in scores.items()) / total_weight
        
    elif method == "geometric_mean":
        # 几何平均(所有得分的乘积开n次方)
        # 避免0值(会导致总体为0)
        non_zero_scores = [max(1e-10, score) for score in scores.values()]
        if len(non_zero_scores) == 0:
            scores["total_score"] = 0.0
            return scores
        
        # 计算几何平均
        composite_score = np.prod(non_zero_scores) ** (1.0 / len(non_zero_scores))
        
    elif method == "harmonic_mean":
        # 调和平均
        # 避免0值(会导致除零错误)
        non_zero_scores = [max(1e-10, score) for score in scores.values()]
        if len(non_zero_scores) == 0:
            scores["total_score"] = 0.0
            return scores
        
        # 计算调和平均
        composite_score = len(non_zero_scores) / sum(1.0 / score for score in non_zero_scores)
        
    elif method == "max":
        # 取最大值
        composite_score = max(scores.values()) if scores else 0.0
        
    elif method == "min":
        # 取最小值
        composite_score = min(scores.values()) if scores else 0.0
        
    else:
        # 默认使用加权平均
        total_weight = sum(weights.get(name, 1.0) for name in scores)
        composite_score = sum(score * weights.get(name, 1.0) for name, score in scores.items()) / total_weight if total_weight > 0 else 0.0
    
    # 确保结果在[0,1]范围内
    composite_score = max(0.0, min(1.0, composite_score))
    
    # 将综合分数添加到结果字典中
    scores["total_score"] = composite_score
    
    return scores


@njit(cache=True, fastmath=True, nogil=True, parallel=True)  
def get_scores_batch(que_list, ref_list=None, tol=(0.003, 0.005)): 
    '''
    (matched_count, matched_ratio, bonanza, simple_dot, modified_dot, entropy) by batch
    ''' 
    n1 = que_list.shape[0]  

    # ref_list为空时，行列都对应que_list  
    if (ref_list is None) or (ref_list.shape[0] == 0):  
        matched_mx = np.zeros((n1, n1), dtype=np.int32)  
        jaccard_mx = np.zeros((n1, n1), dtype=np.float64) 
        bonanza_mx = np.zeros((n1, n1), dtype=np.float64)  
        smp_dot_mx = np.zeros((n1, n1), dtype=np.float64)
        mod_dot_mx = np.zeros((n1, n1), dtype=np.float64)
        entropy_mx = np.zeros((n1, n1), dtype=np.float64)  

        for i in prange(n1):  
            for j in range(i, n1):  
                if i == j:
                    matched_count =  np.sum(que_list[i][:, 0] > 0)  
                    jaccard = bonanza = smp_dot = mod_dot = entropy = 1.0

                else:  
                    matched_count, jaccard, bonanza, smp_dot, mod_dot, entropy = \
                        get_scores_numba(que_list[i], que_list[j], tol) 
                     
                matched_mx[i, j] = matched_mx[j, i] = matched_count 
                jaccard_mx[i, j] = jaccard_mx[j, i] = jaccard 
                bonanza_mx[i, j] = bonanza_mx[j, i] = bonanza  
                smp_dot_mx[i, j] = smp_dot_mx[j, i] = smp_dot 
                mod_dot_mx[i, j] = mod_dot_mx[j, i] = mod_dot
                entropy_mx[i, j] = entropy_mx[j, i] = entropy

    else:  
        n2 = ref_list.shape[0]  
        matched_mx = np.zeros((n1, n2), dtype=np.int32)  
        jaccard_mx = np.zeros((n1, n2), dtype=np.float64) 
        bonanza_mx = np.zeros((n1, n2), dtype=np.float64)  
        smp_dot_mx = np.zeros((n1, n2), dtype=np.float64)
        mod_dot_mx = np.zeros((n1, n2), dtype=np.float64)
        entropy_mx = np.zeros((n1, n2), dtype=np.float64)  

        for i in prange(n1):  
            for j in range(n2):  
                matched_count, jaccard, bonanza, smp_dot, mod_dot, entropy = \
                    get_scores_numba(que_list[i], que_list[j], tol) 
                matched_mx[i, j] = matched_count 
                jaccard_mx[i, j] = jaccard 
                bonanza_mx[i, j] = bonanza  
                smp_dot_mx[i, j] = smp_dot 
                mod_dot_mx[i, j] = mod_dot
                entropy_mx[i, j] = entropy

    return matched_mx, jaccard_mx, bonanza_mx, smp_dot_mx, mod_dot_mx, entropy_mx


@njit(cache=True, fastmath=True, nogil=True)
def get_scores_numba(ms_left, ms_right, tol = (0.003, 0.005)):
    """
    计算两个质谱之间的综合相似度得分 (Numba兼容版本)
    
    参数:
        ms_left: 左谱图数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
        ms_right: 右谱图数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
        tol: 峰匹配的质荷比容差，默认为0.005
        
    返回值:
        综合相似度得分: 范围[0,1]
    """
    # 获取合并峰和质荷比信息
    union_peaks = get_union_ms(ms_left, ms_right, tol[0], tol[1])
    
    # 计算各种相似度得分
    matched_count = get_matched_num(union_peaks)
    matched_ratio = get_matched_peaks_score(union_peaks)
    bonanza = get_bonanza_score(union_peaks)
    simple_dot = get_simple_dot_product(union_peaks)
    modified_dot = get_modified_dot_product_score(union_peaks)
    entropy = get_entropy_similarity(union_peaks)
    return (matched_count, matched_ratio, bonanza, simple_dot, modified_dot, entropy)


@njit(numba.float64[:,:](numba.float64[:,:], numba.float64[:,:], numba.float64, numba.float64),
      cache=True, fastmath=True, nogil=True)
def get_union_ms(ms_left, ms_right, tol1, tol2):
    """
    将两个质谱的峰进行匹配，返回合并后的峰强度数组
    
    参数:
        ms_left: 第一个质谱数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
        ms_right: 第二个质谱数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
        tol1: 母离子容差
        tol2: 碎片离子容差
    
    返回值:
        union_peaks: 形状为(M, 2)的numpy数组，第一列是第一个谱图的峰强度，第二列是第二个谱图的峰强度
    """
    if ms_left.shape[1] != 2 or ms_right.shape[1] != 2:
        raise ValueError("输入数组的形状必须为 (N, 2)")
    
    if ms_left.shape[0] == 0 or ms_right.shape[0] == 0:
        return np.zeros((0, 2), dtype=np.float64)
    
    if abs(ms_left[0, 0] - ms_right[0, 0]) > tol1:
        return np.zeros((0, 2), dtype=np.float64)
    
    left_peaks = ms_left[1:]
    right_peaks = ms_right[1:]
    
    left_peaks_filtered = left_peaks[left_peaks[:, 0] > 0]
    right_peaks_filtered = right_peaks[right_peaks[:, 0] > 0]
    
    if left_peaks_filtered.shape[0] == 0 or right_peaks_filtered.shape[0] == 0:
        return np.zeros((0, 2), dtype=np.float64)
    
    left_peaks_sorted = left_peaks_filtered[np.argsort(left_peaks_filtered[:, 0])]
    right_peaks_sorted = right_peaks_filtered[np.argsort(right_peaks_filtered[:, 0])]
    
    max_size = left_peaks_sorted.shape[0] + right_peaks_sorted.shape[0]
    union_intensity = np.zeros((max_size, 2), dtype=np.float64)
    
    i, j, k = 0, 0, 0
    while i < left_peaks_sorted.shape[0] and j < right_peaks_sorted.shape[0]:
        mz_left = left_peaks_sorted[i, 0]
        mz_right = right_peaks_sorted[j, 0]
        int_left = left_peaks_sorted[i, 1]
        int_right = right_peaks_sorted[j, 1]
        
        if abs(mz_left - mz_right) <= tol2:
            union_intensity[k, 0] = int_left
            union_intensity[k, 1] = int_right
            i += 1
            j += 1
            k += 1
        elif mz_left < mz_right:
            union_intensity[k, 0] = int_left
            union_intensity[k, 1] = 0
            i += 1
            k += 1
        else:
            union_intensity[k, 0] = 0
            union_intensity[k, 1] = int_right
            j += 1
            k += 1
    
    while i < left_peaks_sorted.shape[0]:
        union_intensity[k, 0] = left_peaks_sorted[i, 1]
        union_intensity[k, 1] = 0
        i += 1
        k += 1
    
    while j < right_peaks_sorted.shape[0]:
        union_intensity[k, 0] = 0
        union_intensity[k, 1] = right_peaks_sorted[j, 1]
        j += 1
        k += 1
    
    return union_intensity[:k]


# @njit(numba.float64[:,:](numba.float64[:,:],numba.float64[:,:],numba.float64,numba.float64),
#       cache=True, fastmath=True, nogil=True)
# def get_union_ms(ms_left, ms_right, tol1=0.003, tol2=0.005):
#     """
#     将两个质谱的峰进行匹配，返回合并后的峰强度数组
    
#     参数:
#         ms_left: 第一个质谱数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
#         ms_right: 第二个质谱数据，形状为(N, 2)的numpy数组，第一列是质荷比，第二列是离子丰度
#         tol: 峰匹配的质荷比容差，元组(母离子容差, 碎片离子容差)，默认为(0.003, 0.005)
    
#     返回值:
#         union_peaks: 形状为(M, 2)的numpy数组，第一列是第一个谱图的峰强度，第二列是第二个谱图的峰强度
#     """
#     # 如果任一质谱为空，则返回空数组
#     if ms_left.shape[1] != 2 or ms_right.shape[1] != 2:
#         raise ValueError("The shape of input array must be: (N, 2)")
    
#     # 检查母离子是否匹配
#     if abs(ms_left[0, 0] - ms_right[0, 0]) > tol1:  
#         return np.zeros((0, 2), dtype=np.float64)
            
#     # 跳过母离子，只处理非母离子的峰
#     left_peaks = ms_left[1:].copy()
#     right_peaks = ms_right[1:].copy()
    
#     # 过滤掉质荷比为0的峰
#     left_valid = left_peaks[:, 0] > 0
#     right_valid = right_peaks[:, 0] > 0
    
#     left_peaks_filtered = np.zeros((np.sum(left_valid), 2))
#     right_peaks_filtered = np.zeros((np.sum(right_valid), 2))
    
#     idx = 0
#     for i in range(left_peaks.shape[0]):
#         if left_valid[i]:
#             left_peaks_filtered[idx, 0] = left_peaks[i, 0]
#             left_peaks_filtered[idx, 1] = left_peaks[i, 1]
#             idx += 1
    
#     idx = 0
#     for i in range(right_peaks.shape[0]):
#         if right_valid[i]:
#             right_peaks_filtered[idx, 0] = right_peaks[i, 0]
#             right_peaks_filtered[idx, 1] = right_peaks[i, 1]
#             idx += 1
    
#     left_peaks = left_peaks_filtered
#     right_peaks = right_peaks_filtered
    
#     # 如果过滤后任一质谱为空，则返回空数组
#     if left_peaks.shape[0] == 0 or right_peaks.shape[0] == 0:
#         return np.zeros((0, 2), dtype=np.float64)
    
#     # 按质荷比排序（使用NumPy的argsort和手动重排）
#     left_idx = np.argsort(left_peaks[:, 0])
#     right_idx = np.argsort(right_peaks[:, 0])
    
#     left_peaks_sorted = np.zeros_like(left_peaks)
#     right_peaks_sorted = np.zeros_like(right_peaks)
    
#     for i in range(left_peaks.shape[0]):
#         left_peaks_sorted[i, 0] = left_peaks[left_idx[i], 0]
#         left_peaks_sorted[i, 1] = left_peaks[left_idx[i], 1]
    
#     for i in range(right_peaks.shape[0]):
#         right_peaks_sorted[i, 0] = right_peaks[right_idx[i], 0]
#         right_peaks_sorted[i, 1] = right_peaks[right_idx[i], 1]
    
#     left_peaks = left_peaks_sorted
#     right_peaks = right_peaks_sorted
    
#     # 预分配足够大的数组来存储结果（最坏情况：所有峰都不匹配）
#     max_size = left_peaks.shape[0] + right_peaks.shape[0]
#     union_intensity = np.zeros((max_size, 2))
    
#     # 双指针算法匹配峰
#     i, j, k = 0, 0, 0
#     while i < left_peaks.shape[0] and j < right_peaks.shape[0]:
#         mz_left = left_peaks[i, 0]
#         mz_right = right_peaks[j, 0]
#         int_left = left_peaks[i, 1]
#         int_right = right_peaks[j, 1]
        
#         if abs(mz_left - mz_right) <= tol2:
#             # 匹配峰，添加两个谱图的强度
#             union_intensity[k, 0] = int_left
#             union_intensity[k, 1] = int_right
#             i += 1
#             j += 1
#             k += 1
#         elif mz_left < mz_right:
#             # 左谱图独有峰
#             union_intensity[k, 0] = int_left
#             union_intensity[k, 1] = 0
#             i += 1
#             k += 1
#         else:
#             # 右谱图独有峰
#             union_intensity[k, 0] = 0
#             union_intensity[k, 1] = int_right
#             j += 1
#             k += 1
    
#     # 处理剩余的左谱图峰
#     while i < left_peaks.shape[0]:
#         union_intensity[k, 0] = left_peaks[i, 1]
#         union_intensity[k, 1] = 0
#         i += 1
#         k += 1
    
#     # 处理剩余的右谱图峰
#     while j < right_peaks.shape[0]:
#         union_intensity[k, 0] = 0
#         union_intensity[k, 1] = right_peaks[j, 1]
#         j += 1
#         k += 1
    
#     # 裁剪数组到实际大小
#     return np.ascontiguousarray(union_intensity[:k].astype(np.float64))
