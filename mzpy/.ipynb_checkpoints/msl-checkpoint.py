import numpy as np
from numba import njit

@njit(nogil=True)
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




class MSList: 
    def __init__(self, first_input, second_input=None, precursor_intensity=100):  
        """  
        初始化 MSList_cp 实例，支持两种初始化方式：  
        1. 两个列表: (precursormz_list, msms_list)  
        2. 单一列表: 包含完整谱图数据的列表，或MSList实例  
        
        参数:  
        first_input: 如果提供两个参数，则为母离子的质荷比列表；  
                    如果只提供一个参数，则为完整数据列表或MSList实例  
        second_input: 子离子信息列表，可选  
        precursor_intensity: 母离子的强度，默认为 100  
        """  
        # 检查是否传入了MSList或MSList_cp实例  
        if second_input is None:
            self.data = [np.ascontiguousarray(arr) for arr in first_input]  
        
        # 原始初始化方式：传入两个列表  
        elif second_input is not None:  
            if len(first_input) != len(second_input):  
                raise ValueError("The length of precursormz_list must match the length of msms_list.")
            else:
                self.data = []          
                for precursormz, msms in zip(first_input, second_input):  
                    precursor_array = np.array([[precursormz, precursor_intensity]])         
                    msms = np.asarray(msms)
                    
                    if msms.ndim !=2:
                        raise ValueError(f'Incorrect Data Dimension of msms array: {msms.ndim}')
                        
                    # 将 MS2 数据转为相对丰度后按 mz 从大到小排序
                    max_intensity_ms2 = np.max(msms[:, 1])   
                    if max_intensity_ms2 != 100:  
                        msms[:, 1] = (msms[:, 1] / max_intensity_ms2) * 100  # 转换为相对百分比  
        
                    sorted_msms = msms[np.argsort(msms[:, 0])[::-1]]  # 从大到小排序
        
                    combined_array = np.vstack((precursor_array, sorted_msms))  
                    self.data.append(combined_array)    
        
        else:  
            raise TypeError("Invalid input format. Please provide either: "  
                           "(1) precursormz_list and msms_list, "  
                           "(2) a complete data list, or "  
                           "(3) an MSList instance")
            
        # 不论data的来源如何，self.data最终都转换为numpy数组，保持后续行为的一致性
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
    
    def compute_similarity(self, ref_ms_list, tol=(0.003, 0.005), precursormz_compared=True, n_thread=1):
        if not isinstance(ref_ms_list, self.__class__):
            raise TypeError(f'ref_ms_list must be an instance of {self.__class__.__name__}')

        n_que = len(self.data)
        n_ref = len(ref_ms_list.data)
        similarity_matrix = np.full((n_que, n_ref), -1.0)

        def process_i_chunk(i_start, i_end):
            for i in range(i_start, i_end):
                for j in range(n_ref):
                    similarity_matrix[i, j] = _cosine_(
                        self.data[i], ref_ms_list.data[j], tol, precursormz_compared
                    )

        chunk_size = max(1, (n_que + n_thread - 1) // n_thread)
        chunks = [(start, min(start + chunk_size, n_que)) for start in range(0, n_que, chunk_size)]

        with ThreadPoolExecutor(max_workers=n_thread) as executor:
            futures = []
            for i_start, i_end in chunks:
                futures.append(executor.submit(process_i_chunk, i_start, i_end))
            for future in futures:
                future.result()

        return similarity_matrix

    def compute_similarity_self(self, tol=(0.003, 0.005), precursormz_compared=True, n_thread=1):
        n = len(self.data)
        similarity_matrix = np.full((n, n), -1.0)
        np.fill_diagonal(similarity_matrix, 1.0) # 矩阵对角线

        if n <= 1:
            return similarity_matrix

        def process_i_chunk(i_start, i_end):
            for i in range(i_start, i_end):
                for j in range(i + 1, n):
                    sim = _cosine_(self.data[i], self.data[j], tol, precursormz_compared)
                    similarity_matrix[i, j] = sim
                    similarity_matrix[j, i] = sim

        total_i = n - 1
        chunk_size = max(1, (total_i + n_thread - 1) // n_thread)
        chunks = [(start, min(start + chunk_size, total_i)) for start in range(0, total_i, chunk_size)]

        with ThreadPoolExecutor(max_workers=n_thread) as executor:
            futures = []
            for i_start, i_end in chunks:
                futures.append(executor.submit(process_i_chunk, i_start, i_end))
            for future in futures:
                future.result()

        return similarity_matrix
        

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