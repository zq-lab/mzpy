from typing import List
import numpy as np

class MSdata(np.ndarray):  
    """  
    A custom NumPy array subclass designed for mass spectrometry (MS) data:  
    - Column 0: m/z  
    - Column 1: intensity  
    """  

    def __new__(cls, input_array, metadata=None, to_normalized=True):  
        # Convert the input array to a NumPy array, then view it as MSdata  
        obj = np.asarray(input_array).view(cls)  

        # Verify the shape is (N, 2)  
        if obj.ndim != 2 or obj.shape[1] != 2:  
            raise ValueError("MSdata must be a 2D array with shape (N, 2): [mz, intensity].")  

        # Attach additional metadata if provided  
        obj.metadata = metadata  
        if to_normalized:  
            obj.normalize(inpalce=True)  # Call the normalize method  
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
    
    def filter_out(self, threshold=1):  
        """  
        Filters out rows (after the first) where the intensity (second column) is less than the given threshold.  
        Returns a new MSdata instance with filtered data, keeping the first row.  
        """  
        if self.shape[0] == 0:  
            return self  # Return the empty array if there's no data  
        
        # Keep the first row  
        first_row = self[:1]  
        
        # Create a boolean mask for rows after the first where intensity is >= threshold  
        mask = self[1:, 1] >= threshold  
        
        # Filter the array using the mask and append the first row  
        filtered_data = np.vstack((first_row, self[1:][mask])) if mask.any() else first_row  
        
        # Create a new MSdata object with the filtered data  
        return MSdata(filtered_data, metadata=self.metadata) 
    

    def get_mz(self, intensity=0):
        if intensity == 0:
            return self.mz
        else:
            return self[self[:, 1] > intensity, 0]
        

    @property  
    def mz(self):  
        """  
        Returns the m/z column (column 0 of the array).  
        """  
        return self[:, 0]  
    

    def insert_precursormz(self, mz, intensity=0.0):
        '''
        insert percursor mz at head position
        can not be inserted inpalce
        '''
        return np.insert(self, 0, [mz, intensity], axis=0)
        

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

    def normalize(self, inpalce=False):  
        '''
        Normalize fragment ions to obtain relative intensity and sort them by m/z.
        In-place transformation
        ''' 
        if inpalce:
            if self.shape[0] > 0:
                max_val = self[:, 1].max()  
                if max_val != 0:              # 避免除 0
                    self[:, 1] /= max_val        # 归一化到 0–1，原地修改
                    self[:, 1] *= 100            # 转成百分比，原地修改
        else:
            b = self.copy()
            if b.shape[0] > 0:
                max_val = b[:, 1].max() 
                if max_val != 0:
                    b[:, 1] = b[:, 1] / max_val * 100
            return b


    def to_str(self):
        return str(self.tolist())


