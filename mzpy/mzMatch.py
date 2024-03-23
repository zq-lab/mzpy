import numpy as np
import warnings
from numba import njit
   
def calc_cosine_similarity(a: np.ndarray, b: np.ndarray, tol: float = 0.005) -> float:
    '''
    a: 第一张待求相似度的谱
    b: 第二张待求相似度的谱
    tol: 判定两个峰的质荷比是否相同时的容忍值, 质荷比之差的绝对值小于 tol 时将被视为相等, 默认值为 0.005
    要求输入的谱中各峰按质荷比从小到大排列
    '''
    if not isinstance(a, np.ndarray):
        a = np.array(a)
    if not isinstance(b, np.ndarray):
        b = np.array(b)

    # 创建对齐后的质谱矩阵
    aligned_specs = np.zeros((len(a) + len(b), 2))    # 第 0 列为对齐后的 a 谱强度, 第 1 列为对齐后的 b 谱强度

    # 三个merge指针
    i = 0    # a 中正在被比较的峰的下标
    j = 0    # b 中正在被比较的峰的下标
    k = 0    # aligned_specs 中正在被写入的峰的下标

    while i < len(a) and j < len(b):
        # 如果 a 中正在被比较的峰 小于 b 中正在被比较的峰：
        if a[i, 0] < b[j, 0] - tol:
            # a 列添加该峰，b列对应位置置零
            aligned_specs[k] = (a[i, 1], 0)
            i += 1
        # 如果 a 中正在被比较的峰 大于 b 中正在被比较的峰：
        elif a[i, 0] > b[j, 0] + tol:
            # b 列添加该峰，a列对应位置置零
            aligned_specs[k] = (0, b[j, 1])
            j += 1
        # 如果 a 中正在被比较的峰 约等于 b 中正在被比较的峰：
        else:    # abs(a[i, 0] - b[j, 0]) <= tol
            # a b 列均添加对应峰
            aligned_specs[k] = (a[i, 1], b[j, 1])
            i += 1
            j += 1
        k += 1

    # 如果 a 谱还有未被merge的峰：
    if i < len(a):
        aligned_specs[k: k + len(a) - i, 0] = a[i:, 1]
        k = k + len(a) - i
    # 如果 b 谱还有未被merge的峰：
    elif j < len(b):
        aligned_specs[k: k + len(b) - j, 1] = b[j:, 1]
        k = k + len(b) - j

    # 去掉两列都是0的部分
    aligned_a = aligned_specs[:k, 0]
    aligned_b = aligned_specs[:k, 1]

    # 计算相似度
    similarity = np.dot(aligned_a, aligned_b) / (np.linalg.norm(aligned_a) * np.linalg.norm(aligned_b))

    return similarity