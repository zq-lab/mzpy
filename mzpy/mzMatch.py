import numpy as np
from numba import njit

@njit
def _cosine_(que_sorted, ref_sorted, tol=(0.003, 0.005)):  
    """计算余弦相似度的核心逻辑。
    Assumption: que_sorted and ref_sorted's MS2 are already sorted in descending order by m/z, 
                and relative abundance has been applied.  
    
    Args:  
        que_sorted (ndarray): mz从大到小排序后的查询数据。  
        ref_sorted (ndarray): mz从大到小排序后的参考数据。  
        tol (tuple): 峰匹配的公差，包含两个值。  
    
    Returns:  
        float: 计算得到的余弦相似度或 -1.0。  
    """  
    # 检查 precursor mz 差值  
    mz_diff = abs(que_sorted[0, 0] - ref_sorted[0, 0])  
    if mz_diff > tol[0]:  
        return -1.0  

    # 计算 MS2 的余弦相似度  
    n_que = que_sorted.shape[0]  
    n_ref = ref_sorted.shape[0]  

    union_que = np.empty(n_que + n_ref, dtype=que_sorted.dtype)  
    union_ref = np.empty(n_que + n_ref, dtype=ref_sorted.dtype)  

    i = 1  # 从第二行开始  
    j = 1  # 从第二行开始  
    idx = 0  

    while i < n_que and j < n_ref:  
        mz_que, int_que = que_sorted[i]  
        mz_ref, int_ref = ref_sorted[j]  
        if abs(mz_que - mz_ref) <= tol[1]:  
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
    '''
    union_que and union_ref are two arrays used to store the merged intensity values. 
    Their specific meanings and roles are as follows:

    - union_que: an array used to store the intensity values corresponding to 
                    each matched mz in the query MS (que).
        When a mz in `que` matches a mz in `ref`, `union_que` will store the intensity value in `que`.
    
    The purpose of `union_que` and `union_ref` is to merge the intensity information of the two spectra 
        into a unified structure for subsequent cosine similarity calculation.
    Handling Unmatched Peaks: By setting the unmatched intensity values to 0, 
        it ensures that the cosine similarity calculation is not affected.
    '''
    norm_que = np.sqrt(np.sum(union_que * union_que))  
    norm_ref = np.sqrt(np.sum(union_ref * union_ref))  
    if norm_que == 0 or norm_ref == 0:  
        return 0.0  

    dot_val = np.sum(union_que * union_ref)  
    cos_val = dot_val / (norm_que * norm_ref)  

    lSpectrumCounter = np.sum(union_ref > 10)  

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


class MSmx():
    '''
    Encapsulated MS array, 
        - the first element fixed as the precursor's m/z and intensity (default intensity 100), 
        - the second and subsequent elements as normalized MS2, 
            converted to relative abundance and sorted by m/z in descending order.
    ''' 
    def __init__(self, precursormz, ms2, precursor_intensity=100):  
        """  
        初始化 MSList 实例，创建一个包含母离子和 MS2 数据的数组。  
        
        参数:  
        precursor_intensity: 母离子的强度，默认为 100  
        precursormz: 母离子的质荷比 (float)  
        ms2: MS2 的二维数组，形状为 (N, 2)，表示 [mz, intensity]  
        """  

        # 创建包含母离子质荷比和强度的数组  
        precursor_array = np.array([[precursormz, precursor_intensity]])  

        # 将 MS2 数据转为相对丰度后按 mz 从大到小排序 
        ms2 = np.asarray(ms2) 
        max_intensity_ms2 = np.max(ms2[:, 1])   
        if max_intensity_ms2 != 100:  
            ms2[:, 1] = (ms2[:, 1] / max_intensity_ms2) * 100  # 转换为相对百分比  

        sorted_ms2 = ms2[np.argsort(ms2[:, 0])[::-1]]  # 从大到小排序  

        # 将母离子数组和排序后的 MS2 数据合并  
        data = np.vstack((precursor_array, sorted_ms2))  
        self.data = np.ascontiguousarray(data)


    def __getitem__(self, index):  
        """  
        支持索引访问，返回指定索引的元素。  
        """  
        return self.data[index]  

    def __len__(self):  
        """  
        返回 MSList 中元素的数量。  
        """  
        return len(self.data)  

    def __repr__(self):  
        """  
        返回 MSList 的字符串表示。  
        """  
        return f"MSList({self.data})"

    def compute_similarity(self, ms_mx, tol=(0.003, 0.005)):
        if not isinstance(ms_mx, self.__class__):
            raise TypeError(f'ref_ms_list must be an instance of {self.__class__.__name__}')
        
        if abs(self.precursormz - ms_mx.precursormz) > tol[0]:
            return -1.0
        return _cosine_(self.msms, ms_mx.msms, tol[1])


    @property
    def ms2_mz(self):
        return self.data[1:, 0]

    @property
    def ms2_int(self):
        self.data[1:, 1]

    @property
    def split_ms2(self):
        self.data[1:, 0], self.data[1:, 1]

    @property
    def msms(self):
        return self[1:]    
    
    @property
    def num_ms2(self):
        return self.data.shape[0]-1

    @property
    def precursor(self):
        return self[0]

    @property
    def precursormz(self):
        return self[0][0] 



class MSList:  
    def __init__(self, precursormz_list, msms_list, precursor_intensity=100):
        """  
        初始化 MSList 实例，创建一个包含母离子和对应子离子信息的数组。  
        
        参数:  
        precursormz_list: 母离子的质荷比列表 (list of float)  
        msms_list: 对应的子离子信息列表，每个元素为二维数组 (list of numpy.ndarray)  
        precursor_intensity: 母离子的强度，默认为 100  
        """  
        if len(precursormz_list) != len(msms_list):  
            raise ValueError("The length of precursormz_list must match the length of msms_list.") 

        msms_list = [np.asarray(msms)  for msms in msms_list] 

        # 创建一个空的列表来存储每个母离子及其子离子信息  
        self.data = []  

        for precursormz, msms in zip(precursormz_list, msms_list):  
            # 创建包含母离子质荷比和强度的数组  
            precursor_array = np.array([[precursormz, precursor_intensity]])  
            # 合并母离子数组和子离子信息  

            # 将 MS2 数据转为相对丰度后按 mz 从大到小排序
            max_intensity_ms2 = np.max(msms[:, 1])   
            if max_intensity_ms2 != 100:  
                msms[:, 1] = (msms[:, 1] / max_intensity_ms2) * 100  # 转换为相对百分比  

            sorted_msms = msms[np.argsort(msms[:, 0])[::-1]]  # 从大到小排序

            combined_array = np.vstack((precursor_array, sorted_msms))  
            self.data.append(combined_array)  

        # 将列表转换为 NumPy 数组  
        self.data = np.array(self.data, dtype=object)  # 使用 dtype=object 以支持不同形状的数组  

    def __getitem__(self, index):  
        """  
        支持索引访问，返回指定索引的元素。  
        """  
        return self.data[index]  

    def __len__(self):  
        """  
        返回 MSList 中元素的数量。  
        """  
        return len(self.data)  

    def __repr__(self):  
        """  
        返回 MSList 的字符串表示。  
        """  
        return f"MSList({self.data})" 
    
    def to_cupy(self):  
        """  
        将数据转换为 CuPy 数组并返回。  
        
        返回:  
        list of cupy.ndarray: 每个母离子及其子离子信息的 CuPy 数组列表  
        """  
        import cupy 
        cupy_data = []  
        for ms_array in self.data:  
            # 将 NumPy 数组转换为 CuPy 数组  
            cupy_array = cupy.asarray(ms_array)  # 使用 cp.asarray 进行转换  
            cupy_data.append(cupy_array)  
        return cupy_data 
    
    def to_numpy(self):
        return np.ascontiguousarray(self.data)
    
    def compute_similarity(self, ref_ms_list, tol=(0.003, 0.005)):
        '''
        compoute cosine similarity between self list and ref_ms_list
        ref_ms_list, an instance of MSList
        return:
            a matrix of cosine similarity
        '''
        if not isinstance(ref_ms_list, self.__class__):
            raise TypeError(f'ref_ms_list must be an instance of {self.__class__.__name__}')
        
        n_que = self.data.shape[0]  
        n_ref = ref_ms_list.data.shape[0]  
        similarity_matrix = np.full((n_que, n_ref), -1.0)  # 初始化相似度矩阵，默认值为 -1.0  


        for i in range(n_que):  
            for j in range(n_ref):  
                similarity_matrix[i, j] = self._cosine_(self.data, ref_ms_list.data, tol)  

        return similarity_matrix      
