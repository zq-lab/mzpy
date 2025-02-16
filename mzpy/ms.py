'''
basic transformer or process for MS spectra (ms)

In general ms is a 2D list or numpy matrix (2D array)
    1st column is mz
    2nd column is intensity
'''
import ast
import numpy as np
from typing import List


def to_array(txt, sep1=' ', sep2=':'):
    '''
    transfer ms str to numpy array
    sep1, first level separator
    sep2, second level separator
    '''
    if txt:
        arr = [list(map(float, pair.split(sep2))) for pair in txt.split(sep1)]
        return np.array(arr)
    else:
        return None


def centroid(ms: np.ndarray,
                window_threshold_rate: float=0.33,
                mz_slice_width=0.1,
                n_peaks_threshold = 1) -> List[List[float]]:
    '''
    不同软件的centroid算法结果并不相同
    为了保持一致性，最好使用MS-Dial的centroid结果
    '''
    if len(ms) == 0:
        return []
    if not isinstance(ms, np.ndarray):
        ms = np.array(ms)
    
    uplift = ms[1:] > ms[:-1]
    if not uplift[:, 0].all():
        # 按mz大小排序
        ms = ms[np.argsort(ms[:, 0]), :]
    if len(ms) <= n_peaks_threshold:
        return ms
    
    # 峰检测的向量化操作
    uplift = uplift[:, 1]
    downlift = ms[1:, 1] < ms[:-1, 1]
    peaks_index: List[int] = np.where(uplift[:-1] & downlift[1:])[0] + 1    
    result: List[List[int]] = [None] * peaks_index.shape[0]
    
    for n, pidx in enumerate(peaks_index):
        # 从各峰中心开始，向两侧搜索数据点
        window_size: int = 1                                                        # 搜索的窗口大小
        center_mz, intensity_sum = ms[pidx]                                   # 该峰中心处的 mz, # 该峰中心处的 intensity (用于加权求 mz)
        weighted_mz: float = center_mz * intensity_sum                              # 用于加权求 mz 
        intensity_threshold: float = intensity_sum * window_threshold_rate          # intensity 阈值, 窗口搜索在窗口边界强度低于阈值时结束
        lp: np.ndarray = ms[pidx - 1]     # 窗口左边界的峰
        rp: np.ndarray = ms[pidx + 1]     # 窗口右边界的峰
        
        # 如果:
        # 窗口左边界的峰 intensiy 大于左边界左侧的峰 且
        # 窗口右边界的峰 intensiy 大于右边界右侧的峰 且
        # 窗口左边界与右边界的峰 intensity 均高于 intensity 阈值 且
        # 窗口左边界与右边界的峰 mz 与峰中心 mz 的偏差不超过 mz_slice_width
        # 则向左右扩展窗口        
        while pidx - window_size - 1 >= 0 and \
            pidx + window_size <= peaks_index.shape[0] - 2 and \
            uplift[pidx - window_size - 1] and downlift[pidx + window_size] and \
            (lp := ms[pidx - window_size - 1])[1] > intensity_threshold and \
            (rp := ms[pidx + window_size + 1])[1] > intensity_threshold and \
            abs(lp[0] - center_mz) < mz_slice_width and abs(rp[0] - center_mz) < mz_slice_width:           
            window_size += 1
            intensity_sum += lp[1] + rp[1]
            weighted_mz += lp[0] * lp[1] + rp[0] * rp[1]        
        # 计算加权 mz 后将该峰添加至结果中
        result[n] = [weighted_mz / intensity_sum, ms[pidx][1]]
    
    if not result:
        result: List[List[int]] = [ms[0], ms[-1]]        
    return result


def match_mz(mz1, mz2, tol=0.003, tol_rel=5, mode='abs'):
    '''
    Determine if two mz values (mz1 and mz2) match.
    param:
        mz1, mz2, mz values to compare.
        tol, tolerance. If mode is set as 'both', both tol and tol_rel are required.
    return:
        True or False
    '''
    if mode == 'abs':
        return abs(mz1 - mz2) < tol
    elif mode == 'rel':
        return (1.0E6 * abs(mz1 - mz2) / mz2) < tol
    elif mode == 'both': 
        absolute = abs(mz1-mz2)
        rel = 1.0E6 * absolute / mz2
        return (absolute < tol) or (rel < tol_rel)


def normalize(ms):
    if isinstance(ms, str):
        ms = ast.literal_eval(ms)
    ms = np.ascontiguousarray(ms, dtype=np.float64)
    max_value = ms[:, 1].max()
    if max_value != 100: # 判断是否为相对丰度
        ms[:, 1] = np.around(100 * ms[:, 1] / max_value, decimals=2)        
    ms = ms[ms[:, 0].argsort()] # 按mz排序
    return ms


def to_str(ms):
    if isinstance(ms, np.ndarray):
        ms = ms.tolist()
    return str(ms)