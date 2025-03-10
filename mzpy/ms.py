import ast
from typing import List
import numpy as np
from numba import njit

@njit
def _cosine_(que_sorted, ref_sorted, tol=(0.003, 0.005), precursormz_compared=True):
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
    if precursormz_compared:
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



class MSArray(np.ndarray):  
    """  
    A custom NumPy array subclass designed for mass spectrometry (MS) data:  
    - Column 0: m/z  
    - Column 1: intensity  
    """  

    def __new__(cls, input_array, metadata=None, to_normalized=False):  
        # Convert the input array to a NumPy array, then view it as MSArray 
        #  
        obj = np.asarray(input_array).view(cls)  

        # Verify the shape is (N, 2)  
        if obj.ndim != 2 or obj.shape[1] != 2:  
            raise ValueError("MSArray must be a 2D array with shape (N, 2): [mz, intensity].")  

        # Attach additional metadata if provided  
        obj.metadata = metadata  
        if to_normalized:
            obj = obj.normalize()
        return obj  

    def __array_finalize__(self, obj):  
        """  
        This method is called whenever the array is created (e.g., slicing)  
        to ensure metadata is preserved.  
        """  
        if obj is None:  
            return  
        self.metadata = getattr(obj, 'metadata', None)  


    def centroid(self: np.ndarray,
                    window_threshold_rate: float=0.33,
                    mz_slice_width=0.1,
                    n_peaks_threshold = 1) -> List[List[float]]:
        '''
        不同软件的centroid算法结果并不相同
        为了保持一致性，最好使用MS-Dial的centroid结果
        '''
        if len(self) == 0:
            return []
        if not isinstance(self, np.ndarray):
            self = np.array(self)
        
        uplift = self[1:] > self[:-1]
        if not uplift[:, 0].all():
            # 按mz大小排序
            self = self[np.argsort(self[:, 0]), :]
        if len(self) <= n_peaks_threshold:
            return self
        
        # 峰检测的向量化操作
        uplift = uplift[:, 1]
        downlift = self[1:, 1] < self[:-1, 1]
        peaks_index: List[int] = np.where(uplift[:-1] & downlift[1:])[0] + 1    
        result: List[List[int]] = [None] * peaks_index.shape[0]
        
        for n, pidx in enumerate(peaks_index):
            # 从各峰中心开始，向两侧搜索数据点
            window_size: int = 1                                                        # 搜索的窗口大小
            center_mz, intensity_sum = self[pidx]                                   # 该峰中心处的 mz, # 该峰中心处的 intensity (用于加权求 mz)
            weighted_mz: float = center_mz * intensity_sum                              # 用于加权求 mz 
            intensity_threshold: float = intensity_sum * window_threshold_rate          # intensity 阈值, 窗口搜索在窗口边界强度低于阈值时结束
            lp: np.ndarray = self[pidx - 1]     # 窗口左边界的峰
            rp: np.ndarray = self[pidx + 1]     # 窗口右边界的峰
            
            # 如果:
            # 窗口左边界的峰 intensiy 大于左边界左侧的峰 且
            # 窗口右边界的峰 intensiy 大于右边界右侧的峰 且
            # 窗口左边界与右边界的峰 intensity 均高于 intensity 阈值 且
            # 窗口左边界与右边界的峰 mz 与峰中心 mz 的偏差不超过 mz_slice_width
            # 则向左右扩展窗口        
            while pidx - window_size - 1 >= 0 and \
                pidx + window_size <= peaks_index.shape[0] - 2 and \
                uplift[pidx - window_size - 1] and downlift[pidx + window_size] and \
                (lp := self[pidx - window_size - 1])[1] > intensity_threshold and \
                (rp := self[pidx + window_size + 1])[1] > intensity_threshold and \
                abs(lp[0] - center_mz) < mz_slice_width and abs(rp[0] - center_mz) < mz_slice_width:           
                window_size += 1
                intensity_sum += lp[1] + rp[1]
                weighted_mz += lp[0] * lp[1] + rp[0] * rp[1]        
            # 计算加权 mz 后将该峰添加至结果中
            result[n] = [weighted_mz / intensity_sum, self[pidx][1]]
        
        if not result:
            result: List[List[int]] = [self[0], self[-1]]        
        return result

    @property  
    def mz(self):  
        """  
        Returns the m/z column (column 0 of the array).  
        """  
        return self[:, 0]  

    @property  
    def intensity(self):  
        """  
        Returns the intensity column (column 1 of the array).  
        """  
        return self[:, 1]  

    @property
    def max_intensity_mz(self):  
        """  
        Returns the m/z value corresponding to the maximum intensity.  
        """  
        idx_max = np.argmax(self.intensity)  
        return self.mz[idx_max]

    @property
    def max_mz(self):  
        """  
        Returns the highest m/z value in the current array.  
        """  
        return np.max(self.mz)

    def normalize(self):  
        """  
        Normalizes the intensity to a maximum of 100.  
        Returns a new MSArray instance so that metadata is preserved.  
        """  
        # If the current object is somehow a string, parse it  
        if isinstance(self, str):  
            arr = np.array(ast.literal_eval(self), dtype=np.float64)  
        else:  
            arr = np.ascontiguousarray(self, dtype=np.float64)  

        if arr.shape[0] == 0:  
            # Return an empty MSArray if there's no data  
            empty = np.empty((0,2), dtype=np.float64).view(MSArray)  
            empty.metadata = self.metadata  
            return empty  

        max_value = arr[:, 1].max()  
        if max_value != 100 and max_value != 0:  
            arr[:, 1] = np.around(100.0 * arr[:, 1] / max_value, decimals=2)  

        # Sort by m/z  
        arr = arr[arr[:, 0].argsort()]  

        # Convert back to MSArray to preserve class type and attach metadata  
        normalized_msarray = arr.view(MSArray)  
        normalized_msarray.metadata = self.metadata  
        return normalized_msarray

    def to_str(self):
        return str(self.tolist())



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

        # msms_list = [np.asarray(msms) for msms in msms_list] 

        # 创建一个空的列表来存储每个母离子及其子离子信息  
        self.data = []  

        for precursormz, msms in zip(precursormz_list, msms_list):  
            # 创建包含母离子质荷比和强度的数组  
            precursor_array = np.array([[precursormz, precursor_intensity]])  
            # 合并母离子数组和子离子信息  

            msms = np.asarray(msms)
            # 将 MS2 数据转为相对丰度后按 mz 从大到小排序
            max_intensity_ms2 = np.max(msms[:, 1])   
            if max_intensity_ms2 != 100:  
                msms[:, 1] = (msms[:, 1] / max_intensity_ms2) * 100  # 转换为相对百分比  

            sorted_msms = msms[np.argsort(msms[:, 0])[::-1]]  # 从大到小排序

            combined_array = np.vstack((precursor_array, sorted_msms))  
            self.data.append(combined_array)  

        # 将列表转换为 NumPy 数组  
        self.data = np.ascontiguousarray(self.data, dtype=object)  # 使用 dtype=object 以支持不同形状的数组  

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
       
    def to_numpy(self):
        return np.ascontiguousarray(self.data)
    
    def compute_similarity(self, ref_ms_list, tol=(0.003, 0.005), precursormz_compared=True):
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
                similarity_matrix[i, j] = _cosine_(self.data[i],
                                                   ref_ms_list.data[j],
                                                   tol,
                                                   precursormz_compared)  

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




############    preheat    ############

def preheat_cosine(): 
    print('preheating _cosine_ function') 
    # 创建一些测试数据  
    que_sorted = np.array([[200.0, 100.0], [199.0, 80.0], [198.0, 60.0]])  # 查询数据  
    ref_sorted = np.array([[201.0, 110.0], [200.0, 90.0], [198.5, 70.0]])  # 参考数据  
    
    # 调用该函数以触发JIT编译  
    _ = _cosine_(que_sorted, ref_sorted)  

# 在需要进行计算之前调用预热函数  
preheat_cosine()  
####################################