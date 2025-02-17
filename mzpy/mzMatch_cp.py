import cupy as cp
import numpy as np

from .mzMatch import MSList



def _cosine_(que_sorted, ref_sorted, tol=(0.003, 0.005)):  
    """计算余弦相似度的核心逻辑。  
    Assumption: que_sorted and ref_sorted's MS2 are already sorted in descending order by m/z,   
                and relative abundance has been applied.  
    
    Args:  
        que_sorted (cupy.ndarray): mz从大到小排序后的查询数据。  
        ref_sorted (cupy.ndarray): mz从大到小排序后的参考数据。  
        tol (tuple): 峰匹配的公差，包含两个值。  
    
    Returns:  
        float: 计算得到的余弦相似度或 -1.0。  
    """  
    # 检查 precursor mz 差值  
    mz_diff = cp.abs(que_sorted[0, 0] - ref_sorted[0, 0])  
    if mz_diff > tol[0]:  
        return -1.0  

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



class MSList_cp(MSList):  
    def __init__(self, precursormz_list, msms_list, precursor_intensity=100):
        """  
        初始化 MSList 实例，创建一个包含母离子和对应子离子信息的数组。  
        
        参数:  
        precursormz_list: 母离子的质荷比列表 (list of float)  
        msms_list: 对应的子离子信息列表，每个元素为二维数组 (list of numpy.ndarray)  
        precursor_intensity: 母离子的强度，默认为 100  
        """  
        super().__init__(precursormz_list, msms_list, precursor_intensity)
          

   
        def to_cupy(self):  
            """  
            将数据转换为 CuPy 数组并返回。  
            
            返回:  
            list of cupy.ndarray: 每个母离子及其子离子信息的 CuPy 数组列表  
            """  
            cupy_data = []  
            for ms_array in self.data:  
                # 将 NumPy 数组转换为 CuPy 数组  
                cupy_array = cp.asarray(ms_array)  # 使用 cp.asarray 进行转换  
                cupy_data.append(cupy_array)  
            return cupy_data   
        

    def compute_similarity(self, ref_ms_list, tol=(0.003, 0.005)):  
        '''  
        计算 self list 和 ref_ms_list 之间的余弦相似度  
        ref_ms_list: MSList 的一个实例  
        返回:  
            相似度矩阵  
        '''  
        if not isinstance(ref_ms_list, self.__class__):  
            raise TypeError(f'ref_ms_list must be an instance of {self.__class__.__name__}')  
        
        # 将数据转换为 CuPy 数组  
        self_data = self.to_cupy()  
        ref_data = ref_ms_list.to_cupy()  
        
        n_que = len(self_data)  
        n_ref = len(ref_data)  
        similarity_matrix = cp.full((n_que, n_ref), -1.0)  # 初始化相似度矩阵，默认值为 -1.0  

        for i in range(n_que):  
            for j in range(n_ref):  
                similarity_matrix[i, j] = _cosine_(self_data[i], ref_data[j], tol)  

        return similarity_matrix  
    
    def compute_similarity_self(self, tol=(0.003, 0.005), precursormz_compared=True):
        '''
        self-similarity
        '''       
        n = self.data.shape[0]  
        similarity_matrix = np.full((n, n), -1.0)  # 初始化相似度矩阵，默认值为 -1.0  
        for i in range(n):  
            for j in range(n): 
                if i < j: # 上三角
                    similarity_matrix[i, j] = _cosine_(self.data[i],
                                                       self.data[j],
                                                       tol,
                                                       precursormz_compared) 
                elif i == j: #对角线
                    similarity_matrix[i, j] = 1.0
                elif i > j:
                    similarity_matrix[i, j] = similarity_matrix[j, i]
        return similarity_matrix      


################    preheat    ################         

# 假设这些数据是您准备的母离子和子离子信息  
precursormz_list = [100.0, 150.0, 200.0]  # 示例质荷比列表  
msms_list = [cp.array([[100.0, 50.0], [101.0, 25.0]]),  # 示例子离子信息  
             cp.array([[150.0, 75.0], [151.0, 30.0]]),  
             cp.array([[200.0, 100.0], [201.0, 20.0]])]  

# 创建MSList_cp的实例  
ms_list = MSList_cp(precursormz_list, msms_list)  

# 预热代码  
def preheat_cupy(): 
    print('preheating _cosine_ function')  
    # 创建一个简单的CuPy数组并执行一些操作  
    dummy_array = cp.array([0, 1, 2, 3, 4])  
    _ = cp.sum(dummy_array)  # 执行一个简单的运算以初始化CUDA环境  

# 调用预热函数  
preheat_cupy()  
