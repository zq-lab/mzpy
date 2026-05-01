"""
mzpy.ms
=======

This module provides the core data structure and scoring algorithms for mass
spectrometry (MS and MS/MS) data used in **mzpy**.

Main components
---------------
* **MSdata** – A ``numpy.ndarray`` subclass designed for MS data where
  column 0 holds m/z values and column 1 holds intensities. It supports
  in-place / out-of-place normalization, centroiding, filtering, precursor
  insertion, and stacking of multiple spectra into padded 3-D arrays.

* **Spectral similarity methods** – Python ports of the core scoring
  functions from MS-DIAL (``CompMs.Common.Algorithm.Scoring.MsScanMatching``).
  These methods operate pairwise between the caller spectrum (*self*) and a
  reference spectrum (*ref*):

  - ``simple_dot_product``
  - ``weighted_dot_product``
  - ``reverse_dot_product``
  - ``enhanced_dot_product``
  - ``matched_peaks_scores``
  - ``spectral_entropy_similarity``
  - ``gaussian_similarity``
  - ``isotope_ratio_similarity``

* **Integrated scoring** – Helper functions that combine the individual
  similarity metrics into composite scores exactly as MS-DIAL does:

  - ``total_similarity`` – Multi-dimensional final score combining MS/MS,
    accurate mass, RT, isotope and optional CCS similarities.
  - ``integrated_spectra_similarity`` – Weighted fusion of the three MS/MS
    sub-scores (dot product : reverse dot product : matched peaks = 3:2:1).

All spectral-scoring methods expect both *self* and *ref* to be 2-D arrays
with shape ``(N, 2)`` (m/z, intensity) and assume peaks are sorted by m/z.

* **MSBatch** – A compact CSR-like container for a collection of spectra.
  It can be built directly from a list of ``MSdata`` objects (or nested
  Python lists) or from a zero-padded ``(M, N, 2)`` stack array.
  ``MSBatch`` provides two compute backends selectable at construction:

  - ``device='cpu'`` (default) – exact MS-DIAL dynamic-binning algorithm
    accelerated by Numba (with a pure-Python fallback).
  - ``device='cuda'`` – fixed-grid projection onto a regular m/z lattice.
    Grids are stored as sparse matrices and all algebra runs on the GPU
    via CuPy when available; otherwise it gracefully falls back to CPU
    sparse algebra (still orders of magnitude faster than pairwise loops).

  Matrix methods compute the full M × K similarity matrix in one call:

  - ``simple_dot_product_matrix``
  - ``weighted_dot_product_matrix``
  - ``reverse_dot_product_matrix``
  - ``matched_peaks_scores_matrix``
  - ``enhanced_dot_product_matrix``

  ``prepare_grid(bin, mass_begin, mass_end)`` exposes the internal sparse
  grid representation for custom downstream computations.
"""

from typing import List, Sequence, Union
import numpy as np

try:
    import numba as nb
    _HAS_NUMBA = True
except Exception:
    _HAS_NUMBA = False

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


    @staticmethod
    def stack(ms_list, n=None):
        """
        Stack multiple MSdata or 2-D arrays into a padded 3-D array of shape (M, N, 2).

        Each element in *ms_list* should have shape (K, 2).  The resulting array
        has shape ``(len(ms_list), N, 2)`` where ``N = max(K)`` unless explicitly
        provided by *n*.  Shorter arrays are zero-padded at the end.

        Parameters
        ----------
        ms_list : list of array-like
            List of MSdata instances or (K, 2) arrays.
        n : int, optional
            Fixed number of rows N.  If None, uses the maximum length found in
            *ms_list*.

        Returns
        -------
        np.ndarray
            Padded 3-D array with dtype float64 and shape (M, N, 2).
        """
        if not ms_list:
            return np.empty((0, 0, 2), dtype=np.float64)

        arrays = [np.asarray(ms) for ms in ms_list]

        if n is None:
            n = max(arr.shape[0] for arr in arrays)

        m = len(arrays)
        result = np.zeros((m, n, 2), dtype=np.float64)

        for i, arr in enumerate(arrays):
            k = min(arr.shape[0], n)
            if k > 0:
                result[i, :k, :] = arr[:k, :]

        return result


    def to_str(self):
        return str(self.tolist())

    # ------------------------------------------------------------------
    #  MS-DIAL spectral similarity/scoring methods (ported from
    #  CompMs.Common.Algorithm.Scoring.MsScanMatching)
    # ------------------------------------------------------------------

    def _get_sorted_copy(self):
        """Return a copy sorted by m/z if the current array is not sorted."""
        if self.shape[0] > 1 and np.any(self[1:, 0] < self[:-1, 0]):
            return self[np.argsort(self[:, 0]), :]
        return self

    def simple_dot_product(self, ref, bin: float = 0.05,
                           mass_begin: float = 0.0, mass_end: float = 2000.0) -> float:
        """
        Simple dot-product similarity (squared) between *self* and *ref*.
        Peaks are binned by ``bin`` (Da) and normalized to the base-peak
        intensity before the cosine-like calculation.

        Reference
        ---------
        Stein, S. E. *J. Am. Soc. Mass. Spectrom.* 1999, 10, 770–781.
        """
        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return -1.0

        peaks1 = peaks1[peaks1[:, 0].argsort()] if peaks1.shape[0] > 1 else peaks1
        peaks2 = peaks2[peaks2[:, 0].argsort()] if peaks2.shape[0] > 1 else peaks2

        max_mz = min(mass_end, max(peaks1[-1, 0], peaks2[-1, 0]))
        remaind_m, remaind_l = 0, 0

        while remaind_m < peaks1.shape[0] and peaks1[remaind_m, 0] < mass_begin - bin:
            remaind_m += 1
        while remaind_l < peaks2.shape[0] and peaks2[remaind_l, 0] < mass_begin - bin:
            remaind_l += 1

        if remaind_m >= peaks1.shape[0] or remaind_l >= peaks2.shape[0]:
            return 0.0

        focused_mz = min(peaks1[remaind_m, 0], peaks2[remaind_l, 0])
        measured, reference = [], []
        base_m = base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < peaks1.shape[0]:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < peaks2.shape[0]:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            measured.append([focused_mz, sum_m])
            reference.append([focused_mz, sum_r])
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin > max(peaks1[-1, 0], peaks2[-1, 0]):
                break
            if remaind_m >= peaks1.shape[0] or remaind_l >= peaks2.shape[0]:
                focused_mz = (peaks1[remaind_m, 0] if remaind_l >= peaks2.shape[0]
                              else peaks2[remaind_l, 0])
                continue
            next_m = peaks1[remaind_m, 0]
            next_l = peaks2[remaind_l, 0]
            if focused_mz + bin > next_l and focused_mz + bin <= next_m:
                focused_mz = next_m
            elif focused_mz + bin <= next_l and focused_mz + bin > next_m:
                focused_mz = next_l
            else:
                focused_mz = min(next_m, next_l)

        if base_m == 0 or base_r == 0:
            return 0.0

        scalar_m = scalar_r = covariance = 0.0
        for i in range(len(measured)):
            m_int = measured[i][1] / base_m * 999.0
            r_int = reference[i][1] / base_r * 999.0
            scalar_m += m_int
            scalar_r += r_int
            covariance += np.sqrt(m_int * r_int)

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r

    def weighted_dot_product(self, ref, bin: float = 0.05,
                             mass_begin: float = 0.0, mass_end: float = 2000.0) -> float:
        """
        Weighted dot-product similarity (squared) between *self* and *ref*.
        Unlike the simple dot product, intensities are weighted by m/z and a
        peak-count penalty is applied for very small reference spectra.
        """
        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return -1.0

        peaks1 = peaks1[peaks1[:, 0].argsort()] if peaks1.shape[0] > 1 else peaks1
        peaks2 = peaks2[peaks2[:, 0].argsort()] if peaks2.shape[0] > 1 else peaks2

        min_mz = max(mass_begin, min(peaks1[0, 0], peaks2[0, 0]))
        max_mz = min(mass_end, max(peaks1[-1, 0], peaks2[-1, 0]))
        focused_mz = min_mz
        remaind_m = remaind_l = 0

        measured, reference = [], []
        base_m = base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < peaks1.shape[0]:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < peaks2.shape[0]:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            measured.append([focused_mz, sum_m])
            reference.append([focused_mz, sum_r])
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin > max(peaks1[-1, 0], peaks2[-1, 0]):
                break
            next_m = peaks1[remaind_m, 0] if remaind_m < peaks1.shape[0] else float('inf')
            next_l = peaks2[remaind_l, 0] if remaind_l < peaks2.shape[0] else float('inf')
            if focused_mz + bin > next_l and focused_mz + bin <= next_m:
                focused_mz = next_m
            elif focused_mz + bin <= next_l and focused_mz + bin > next_m:
                focused_mz = next_l
            else:
                focused_mz = min(next_m, next_l)

        if base_m == 0 or base_r == 0:
            return 0.0

        sum_measure = sum_reference = 0.0
        e_counter = l_counter = 0
        for i in range(len(measured)):
            m_norm = measured[i][1] / base_m
            r_norm = reference[i][1] / base_r
            sum_measure += m_norm
            sum_reference += r_norm
            if m_norm > 0.1:
                e_counter += 1
            if r_norm > 0.1:
                l_counter += 1

        peak_count_penalty = 1.0
        if l_counter == 1:
            peak_count_penalty = 0.75
        elif l_counter == 2:
            peak_count_penalty = 0.88
        elif l_counter == 3:
            peak_count_penalty = 0.94
        elif l_counter == 4:
            peak_count_penalty = 0.97

        scalar_m = scalar_r = covariance = 0.0
        cutoff = 0.01
        for i in range(len(measured)):
            m_int = measured[i][1] / base_m
            if m_int < cutoff:
                continue
            r_int = reference[i][1] / base_r
            mz = measured[i][0]
            scalar_m += m_int * mz
            scalar_r += r_int * mz
            covariance += np.sqrt(m_int * r_int) * mz

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r * peak_count_penalty

    def reverse_dot_product(self, ref, bin: float = 0.05,
                            mass_begin: float = 0.0, mass_end: float = 2000.0) -> float:
        """
        Reverse dot-product similarity (squared) between *self* and *ref*.
        Only bins that contain peaks in the *reference* spectrum are scored;
        unmatched query peaks are ignored.
        """
        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return -1.0

        peaks1 = peaks1[peaks1[:, 0].argsort()] if peaks1.shape[0] > 1 else peaks1
        peaks2 = peaks2[peaks2[:, 0].argsort()] if peaks2.shape[0] > 1 else peaks2

        min_mz = max(mass_begin, peaks2[0, 0])
        max_mz = min(mass_end, peaks2[-1, 0])
        focused_mz = min_mz
        remaind_m = remaind_l = 0

        measured, reference = [], []
        base_m = base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < peaks1.shape[0]:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < peaks2.shape[0]:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            measured.append([focused_mz, sum_m])
            reference.append([focused_mz, sum_r])
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin > peaks2[-1, 0]:
                break
            if remaind_l < peaks2.shape[0]:
                focused_mz = peaks2[remaind_l, 0]
            else:
                break

        if base_m == 0 or base_r == 0:
            return 0.0

        sum_measure = sum_reference = 0.0
        e_counter = l_counter = 0
        for i in range(len(measured)):
            m_norm = measured[i][1] / base_m
            r_norm = reference[i][1] / base_r
            sum_measure += m_norm
            sum_reference += r_norm
            if m_norm > 0.1:
                e_counter += 1
            if r_norm > 0.1:
                l_counter += 1

        peak_count_penalty = 1.0
        if l_counter == 1:
            peak_count_penalty = 0.75
        elif l_counter == 2:
            peak_count_penalty = 0.88
        elif l_counter == 3:
            peak_count_penalty = 0.94
        elif l_counter == 4:
            peak_count_penalty = 0.97

        scalar_m = scalar_r = covariance = 0.0
        cutoff = 0.01
        for i in range(len(measured)):
            r_int = reference[i][1] / base_r
            if r_int < cutoff:
                continue
            m_int = measured[i][1] / base_m
            mz = measured[i][0]
            scalar_m += m_int * mz
            scalar_r += r_int * mz
            covariance += np.sqrt(m_int * r_int) * mz

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r * peak_count_penalty

    def matched_peaks_scores(self, ref, bin: float = 0.05,
                             mass_begin: float = 0.0, mass_end: float = 2000.0):
        """
        Matched-peaks percentage and count.

        Returns
        -------
        list
            ``[matched_ratio, matched_count]`` where *matched_ratio* is the
            fraction of reference peaks (above 1 % of base-peak) found in the
            query spectrum.
        """
        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return [-1, -1]

        peaks1 = peaks1[peaks1[:, 0].argsort()] if peaks1.shape[0] > 1 else peaks1
        peaks2 = peaks2[peaks2[:, 0].argsort()] if peaks2.shape[0] > 1 else peaks2

        min_mz = max(peaks2[0, 0], mass_begin)
        max_mz = min(peaks2[-1, 0], mass_end)
        focused_mz = min_mz
        remaind_m = remaind_l = 0
        max_lib_intensity = float(np.max(peaks2[:, 1]))
        counter = 0
        lib_counter = 0

        while focused_mz <= max_mz:
            sum_l = 0.0
            i = remaind_l
            while i < peaks2.shape[0]:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_l += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            if sum_l >= 0.01 * max_lib_intensity:
                lib_counter += 1

            sum_m = 0.0
            i = remaind_m
            while i < peaks1.shape[0]:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            if sum_m > 0 and sum_l >= 0.01 * max_lib_intensity:
                counter += 1

            if focused_mz + bin > peaks2[-1, 0]:
                break
            if remaind_l < peaks2.shape[0]:
                focused_mz = peaks2[remaind_l, 0]
            else:
                break

        if lib_counter == 0:
            return [0.0, 0]
        return [counter / lib_counter, counter]

    def enhanced_dot_product(self, ref, bin: float = 0.05,
                             mass_begin: float = 0.0, mass_end: float = 2000.0,
                             penalty: float = 0.5) -> float:
        """
        Enhanced dot-product similarity (squared).

        Parameters
        ----------
        penalty : float, default 0.5
            Weight for unmatched query peaks.  ``penalty=1`` behaves like the
            reverse dot product; ``penalty=0`` behaves like the standard dot
            product.

        Reference
        ---------
        Xing, S. et al. *Anal. Chem.* 2025, 97(33).
        """
        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return -1.0

        peaks1 = peaks1[peaks1[:, 0].argsort()] if peaks1.shape[0] > 1 else peaks1
        peaks2 = peaks2[peaks2[:, 0].argsort()] if peaks2.shape[0] > 1 else peaks2

        min_mz = max(mass_begin, min(peaks1[0, 0], peaks2[0, 0]))
        max_mz = min(mass_end, max(peaks1[-1, 0], peaks2[-1, 0]))
        focused_mz = min_mz
        remaind_m = remaind_l = 0

        measured, reference = [], []
        base_m = base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < peaks1.shape[0]:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < peaks2.shape[0]:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin:
                    i += 1
                    continue
                if focused_mz - bin <= mz < focused_mz + bin:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            measured.append([focused_mz, sum_m])
            reference.append([focused_mz, sum_r])
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin > max(peaks1[-1, 0], peaks2[-1, 0]):
                break
            next_m = peaks1[remaind_m, 0] if remaind_m < peaks1.shape[0] else float('inf')
            next_l = peaks2[remaind_l, 0] if remaind_l < peaks2.shape[0] else float('inf')
            if focused_mz + bin > next_l and focused_mz + bin <= next_m:
                focused_mz = next_m
            elif focused_mz + bin <= next_l and focused_mz + bin > next_m:
                focused_mz = next_l
            else:
                focused_mz = min(next_m, next_l)

        if base_m == 0 or base_r == 0:
            return 0.0

        e_counter = l_counter = 0
        for i in range(len(measured)):
            m_norm = measured[i][1] / base_m
            r_norm = reference[i][1] / base_r
            if m_norm > 0.1:
                e_counter += 1
            if r_norm > 0.1:
                l_counter += 1

        peak_count_penalty = 1.0
        if l_counter == 1:
            peak_count_penalty = 0.75
        elif l_counter == 2:
            peak_count_penalty = 0.88
        elif l_counter == 3:
            peak_count_penalty = 0.94
        elif l_counter == 4:
            peak_count_penalty = 0.97

        scalar_m = scalar_r = covariance = 0.0
        cutoff = 0.01
        for i in range(len(measured)):
            r_int = reference[i][1] / base_r
            if r_int < cutoff:
                continue
            m_int = measured[i][1] / base_m
            mz = measured[i][0]
            scalar_r += r_int * mz
            if m_int == 0.0:
                scalar_m += m_int * (1.0 - penalty) * mz
            else:
                scalar_m += m_int * mz
            covariance += np.sqrt(m_int * r_int) * mz

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r * peak_count_penalty

    def spectral_entropy_similarity(self, ref, bin: float = 0.05) -> float:
        """
        Spectral-entropy similarity between *self* and *ref*.

        Returns
        -------
        float
            1 - (2*H12 - H1 - H2) / 2,  where H is the Shannon entropy of the
            intensity distribution.
        """
        def _norm_total(p):
            s = p[:, 1].sum()
            if s == 0:
                return p.copy()
            out = p.copy()
            out[:, 1] /= s
            return out

        def _binned(p, b):
            if p.shape[0] == 0:
                return p.copy()
            p = p[p[:, 0].argsort()]
            res = []
            i = 0
            while i < p.shape[0]:
                mz = p[i, 0]
                intens = 0.0
                while i < p.shape[0] and abs(p[i, 0] - mz) < b:
                    intens += p[i, 1]
                    i += 1
                res.append([mz, intens])
            return np.array(res)

        def _combined(p1, p2, b):
            if p1.shape[0] == 0 and p2.shape[0] == 0:
                return np.empty((0, 2))
            c = np.vstack([p1, p2])
            c = c[c[:, 0].argsort()]
            res = []
            i = 0
            while i < c.shape[0]:
                mz = c[i, 0]
                intens = 0.0
                while i < c.shape[0] and abs(c[i, 0] - mz) < b:
                    intens += c[i, 1]
                    i += 1
                res.append([mz, intens])
            return np.array(res)

        def _entropy(p):
            s = p[:, 1].sum()
            if s == 0:
                return 0.0
            prob = p[:, 1] / s
            prob = prob[prob > 0]
            return float(-np.sum(prob * np.log2(prob)))

        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return -1.0

        norm1 = _norm_total(peaks1)
        norm2 = _norm_total(peaks2)
        binned1 = _binned(norm1, bin)
        binned2 = _binned(norm2, bin)
        combined = _combined(norm1, norm2, bin)

        h12 = _entropy(combined)
        h1 = _entropy(binned1)
        h2 = _entropy(binned2)
        return 1.0 - (2.0 * h12 - h1 - h2) * 0.5

    def gaussian_similarity(self, ref, tolerance: float = 0.01) -> float:
        """
        Gaussian similarity for a single scalar value (e.g. precursor m/z).
        Compares the first m/z value of *self* and *ref*.

        Reference
        ---------
        Tsugawa, H. et al. *Anal. Chem.* 2013, 85, 5191–5199.
        """
        actual = float(self[0, 0]) if self.shape[0] > 0 else 0.0
        reference = float(ref[0, 0]) if np.asarray(ref).shape[0] > 0 else 0.0
        if actual <= 0 or reference <= 0:
            return -1.0
        return float(np.exp(-0.5 * ((actual - reference) / tolerance) ** 2))

    def isotope_ratio_similarity(self, ref, targeted_mz: float = None,
                                  tolerance: float = 0.01) -> float:
        """
        Isotope-ratio similarity between *self* and *ref*.
        Both arrays should contain isotopic peaks (m/z, relative abundance).

        Parameters
        ----------
        targeted_mz : float, optional
            Precursor m/z of the query.  Defaults to ``self[0, 0]``.
        """
        peaks1 = np.asarray(self)
        peaks2 = np.asarray(ref)
        if peaks1.shape[0] == 0 or peaks2.shape[0] == 0:
            return -1.0

        if peaks1[0, 1] <= 0 or peaks2[0, 1] <= 0:
            return -1.0

        n = min(peaks1.shape[0], peaks2.shape[0])
        similarity = 0.0
        for i in range(1, n):
            ratio1 = peaks1[i, 1] / peaks1[0, 1]
            ratio2 = peaks2[i, 1] / peaks2[0, 1]
            if ratio1 <= 1 and ratio2 <= 1:
                similarity += abs(ratio1 - ratio2)
            elif ratio1 > ratio2:
                similarity += 1.0 - ratio2 / ratio1
            else:
                similarity += 1.0 - ratio1 / ratio2
        return 1.0 - similarity

    @staticmethod
    def total_similarity(accurate_mass_similarity: float,
                         rt_similarity: float,
                         isotope_similarity: float,
                         spectra_similarity: float,
                         reverse_search_similarity: float,
                         presence_similarity: float,
                         spectrum_penalty: bool = False,
                         target_omics: str = 'metabolomics',
                         is_use_rt: bool = True,
                         ccs_similarity: float = None) -> float:
        """
        MS-DIAL ``GetTotalSimilarity`` 的 Python 移植版本。

        将 MS/MS 光谱相似度、精确质量、保留时间、同位素（以及可选的 CCS）
        等多个维度的分数按 MS-DIAL 内置权重融合为单一的最终排名分数。

        Parameters
        ----------
        accurate_mass_similarity : float
            精确质量匹配分数（如 ``gaussian_similarity`` 的结果）。
        rt_similarity : float
            保留时间相似度（< 0 表示未使用或不可用）。
        isotope_similarity : float
            同位素比例相似度（< 0 表示未使用或不可用）。
        spectra_similarity : float
            光谱相似度（通常取 ``weighted_dot_product`` 的平方根）。
        reverse_search_similarity : float
            反向搜索相似度（通常取 ``reverse_dot_product`` 的平方根）。
        presence_similarity : float
            匹配峰百分比（``matched_peaks_scores`` 返回的 ratio）。
        spectrum_penalty : bool, default False
            若为 True 且 target_omics='metabolomics'，MS/MS 分数会乘以 0.5 的惩罚。
        target_omics : {'metabolomics', 'lipidomics'}, default 'metabolomics'
            目标组学类型，决定各维度的权重。
        is_use_rt : bool, default True
            是否在总分中纳入保留时间分数。
        ccs_similarity : float, optional
            CCS 相似度。若提供，则启用 CCS 维度并按含 CCS 的公式计算。

        Returns
        -------
        float
            综合总相似度分数（0 ~ 1 范围，但可能因输入值而超出）。
        """
        omics = target_omics.lower()
        if omics == 'lipidomics':
            dot_product_factor = 1.0
            reverse_dot_prod_factor = 2.0
            presence_percentage_factor = 3.0
            msms_factor = 1.5
            rt_factor = 0.5
            mass_factor = 1.0
            isotope_factor = 0.0
            ccs_factor = 1.0
        else:  # metabolomics
            dot_product_factor = 3.0
            reverse_dot_prod_factor = 2.0
            presence_percentage_factor = 1.0
            msms_factor = 2.0
            rt_factor = 1.0
            mass_factor = 1.0
            isotope_factor = 0.0
            ccs_factor = 2.0

        msms_similarity = (
            dot_product_factor * spectra_similarity
            + reverse_dot_prod_factor * reverse_search_similarity
            + presence_percentage_factor * presence_similarity
        ) / (dot_product_factor + reverse_dot_prod_factor + presence_percentage_factor)

        if spectrum_penalty and omics == 'metabolomics':
            msms_similarity *= 0.5

        # 使用含 CCS 的分支逻辑（与 C# 含 CCS 重载保持一致）
        if ccs_similarity is not None:
            use_rt = is_use_rt and rt_similarity >= 0
            use_ccs = ccs_similarity >= 0
            use_isotope = isotope_similarity >= 0

            numerator = msms_factor * msms_similarity + mass_factor * accurate_mass_similarity
            denominator = msms_factor + mass_factor

            if use_rt:
                numerator += rt_factor * rt_similarity
                denominator += rt_factor
            if use_ccs:
                numerator += ccs_factor * ccs_similarity
                denominator += ccs_factor
            if use_isotope:
                numerator += isotope_factor * isotope_similarity
                denominator += isotope_factor

            return numerator / denominator

        # 不含 CCS 的分支逻辑（与 C# 无 CCS 重载保持一致）
        if not is_use_rt:
            if isotope_similarity < 0:
                return (msms_factor * msms_similarity + mass_factor * accurate_mass_similarity) \
                       / (msms_factor + mass_factor)
            else:
                return (msms_factor * msms_similarity + mass_factor * accurate_mass_similarity
                        + isotope_factor * isotope_similarity) \
                       / (msms_factor + mass_factor + isotope_factor)
        else:
            if rt_similarity < 0 and isotope_similarity < 0:
                return (msms_factor * msms_similarity + mass_factor * accurate_mass_similarity) \
                       / (msms_factor + mass_factor + rt_factor)
            elif rt_similarity < 0 and isotope_similarity >= 0:
                return (msms_factor * msms_similarity + mass_factor * accurate_mass_similarity
                        + isotope_factor * isotope_similarity) \
                       / (msms_factor + mass_factor + isotope_factor + rt_factor)
            elif isotope_similarity < 0 and rt_similarity >= 0:
                return (msms_factor * msms_similarity + mass_factor * accurate_mass_similarity
                        + rt_factor * rt_similarity) \
                       / (msms_factor + mass_factor + rt_factor)
            else:
                return (msms_factor * msms_similarity + mass_factor * accurate_mass_similarity
                        + rt_factor * rt_similarity + isotope_factor * isotope_similarity) \
                       / (msms_factor + mass_factor + rt_factor + isotope_factor)

    @staticmethod
    def integrated_spectra_similarity(spectra_similarity: float,
                                       reverse_search_similarity: float,
                                       presence_similarity: float) -> float:
        """
        MS-DIAL ``GetIntegratedSpectraSimilarity`` 的 Python 移植版本。

        仅对 MS/MS 层面的三个子分数进行加权综合，不涉及 mass/RT/isotope 等维度。

        Parameters
        ----------
        spectra_similarity : float
            光谱相似度（如 weighted dot product 的平方根）。
        reverse_search_similarity : float
            反向搜索相似度（如 reverse dot product 的平方根）。
        presence_similarity : float
            匹配峰百分比。

        Returns
        -------
        float
            综合光谱相似度。
        """
        dotproduct_fact = 3.0
        rev_dotproduct_fact = 2.0
        matched_ratio_fact = 1.0
        return (dotproduct_fact * spectra_similarity
                + rev_dotproduct_fact * reverse_search_similarity
                + matched_ratio_fact * presence_similarity) \
               / (dotproduct_fact + rev_dotproduct_fact + matched_ratio_fact)


# ---------------------------------------------------------------------------
#  Numba kernels & pure-Python fallbacks for batch spectral similarity
#  (ported from MS-DIAL MsScanMatching)
# ---------------------------------------------------------------------------

if _HAS_NUMBA:
    @nb.njit
    def _simple_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        n1 = peaks1.shape[0]
        n2 = peaks2.shape[0]
        if n1 == 0 or n2 == 0:
            return 0.0

        max_mz = min(mass_end, max(peaks1[n1 - 1, 0], peaks2[n2 - 1, 0]))
        remaind_m = 0
        remaind_l = 0

        while remaind_m < n1 and peaks1[remaind_m, 0] < mass_begin - bin_width:
            remaind_m += 1
        while remaind_l < n2 and peaks2[remaind_l, 0] < mass_begin - bin_width:
            remaind_l += 1

        if remaind_m >= n1 or remaind_l >= n2:
            return 0.0

        focused_mz = min(peaks1[remaind_m, 0], peaks2[remaind_l, 0])
        max_size = n1 + n2
        sum_m_arr = np.empty(max_size, dtype=np.float64)
        sum_r_arr = np.empty(max_size, dtype=np.float64)
        size = 0
        base_m = 0.0
        base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < n1:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < n2:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            sum_m_arr[size] = sum_m
            sum_r_arr[size] = sum_r
            size += 1
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin_width > max(peaks1[n1 - 1, 0], peaks2[n2 - 1, 0]):
                break
            if remaind_m >= n1 or remaind_l >= n2:
                focused_mz = (peaks1[remaind_m, 0] if remaind_l >= n2
                              else peaks2[remaind_l, 0])
                continue
            next_m = peaks1[remaind_m, 0]
            next_l = peaks2[remaind_l, 0]
            if focused_mz + bin_width > next_l and focused_mz + bin_width <= next_m:
                focused_mz = next_m
            elif focused_mz + bin_width <= next_l and focused_mz + bin_width > next_m:
                focused_mz = next_l
            else:
                focused_mz = min(next_m, next_l)

        if base_m == 0 or base_r == 0:
            return 0.0

        scalar_m = 0.0
        scalar_r = 0.0
        covariance = 0.0
        for i in range(size):
            m_int = sum_m_arr[i] / base_m * 999.0
            r_int = sum_r_arr[i] / base_r * 999.0
            scalar_m += m_int
            scalar_r += r_int
            covariance += np.sqrt(m_int * r_int)

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r

    @nb.njit
    def _weighted_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        n1 = peaks1.shape[0]
        n2 = peaks2.shape[0]
        if n1 == 0 or n2 == 0:
            return 0.0

        min_mz = max(mass_begin, min(peaks1[0, 0], peaks2[0, 0]))
        max_mz = min(mass_end, max(peaks1[n1 - 1, 0], peaks2[n2 - 1, 0]))
        focused_mz = min_mz
        remaind_m = 0
        remaind_l = 0

        max_size = n1 + n2
        mz_arr = np.empty(max_size, dtype=np.float64)
        sum_m_arr = np.empty(max_size, dtype=np.float64)
        sum_r_arr = np.empty(max_size, dtype=np.float64)
        size = 0
        base_m = 0.0
        base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < n1:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < n2:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            mz_arr[size] = focused_mz
            sum_m_arr[size] = sum_m
            sum_r_arr[size] = sum_r
            size += 1
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin_width > max(peaks1[n1 - 1, 0], peaks2[n2 - 1, 0]):
                break
            next_m = peaks1[remaind_m, 0] if remaind_m < n1 else np.inf
            next_l = peaks2[remaind_l, 0] if remaind_l < n2 else np.inf
            if focused_mz + bin_width > next_l and focused_mz + bin_width <= next_m:
                focused_mz = next_m
            elif focused_mz + bin_width <= next_l and focused_mz + bin_width > next_m:
                focused_mz = next_l
            else:
                focused_mz = min(next_m, next_l)

        if base_m == 0 or base_r == 0:
            return 0.0

        sum_measure = 0.0
        sum_reference = 0.0
        l_counter = 0
        for i in range(size):
            m_norm = sum_m_arr[i] / base_m
            r_norm = sum_r_arr[i] / base_r
            sum_measure += m_norm
            sum_reference += r_norm
            if r_norm > 0.1:
                l_counter += 1

        peak_count_penalty = 1.0
        if l_counter == 1:
            peak_count_penalty = 0.75
        elif l_counter == 2:
            peak_count_penalty = 0.88
        elif l_counter == 3:
            peak_count_penalty = 0.94
        elif l_counter == 4:
            peak_count_penalty = 0.97

        scalar_m = 0.0
        scalar_r = 0.0
        covariance = 0.0
        cutoff = 0.01
        for i in range(size):
            m_int = sum_m_arr[i] / base_m
            if m_int < cutoff:
                continue
            r_int = sum_r_arr[i] / base_r
            mz = mz_arr[i]
            scalar_m += m_int * mz
            scalar_r += r_int * mz
            covariance += np.sqrt(m_int * r_int) * mz

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r * peak_count_penalty

    @nb.njit
    def _reverse_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        n1 = peaks1.shape[0]
        n2 = peaks2.shape[0]
        if n1 == 0 or n2 == 0:
            return 0.0

        min_mz = max(mass_begin, peaks2[0, 0])
        max_mz = min(mass_end, peaks2[n2 - 1, 0])
        focused_mz = min_mz
        remaind_m = 0
        remaind_l = 0

        max_size = n1 + n2
        mz_arr = np.empty(max_size, dtype=np.float64)
        sum_m_arr = np.empty(max_size, dtype=np.float64)
        sum_r_arr = np.empty(max_size, dtype=np.float64)
        size = 0
        base_m = 0.0
        base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < n1:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < n2:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            mz_arr[size] = focused_mz
            sum_m_arr[size] = sum_m
            sum_r_arr[size] = sum_r
            size += 1
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin_width > peaks2[n2 - 1, 0]:
                break
            if remaind_l < n2:
                focused_mz = peaks2[remaind_l, 0]
            else:
                break

        if base_m == 0 or base_r == 0:
            return 0.0

        sum_measure = 0.0
        sum_reference = 0.0
        l_counter = 0
        for i in range(size):
            m_norm = sum_m_arr[i] / base_m
            r_norm = sum_r_arr[i] / base_r
            sum_measure += m_norm
            sum_reference += r_norm
            if r_norm > 0.1:
                l_counter += 1

        peak_count_penalty = 1.0
        if l_counter == 1:
            peak_count_penalty = 0.75
        elif l_counter == 2:
            peak_count_penalty = 0.88
        elif l_counter == 3:
            peak_count_penalty = 0.94
        elif l_counter == 4:
            peak_count_penalty = 0.97

        scalar_m = 0.0
        scalar_r = 0.0
        covariance = 0.0
        cutoff = 0.01
        for i in range(size):
            r_int = sum_r_arr[i] / base_r
            if r_int < cutoff:
                continue
            m_int = sum_m_arr[i] / base_m
            mz = mz_arr[i]
            scalar_m += m_int * mz
            scalar_r += r_int * mz
            covariance += np.sqrt(m_int * r_int) * mz

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r * peak_count_penalty

    @nb.njit
    def _matched_peaks_scores_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        n1 = peaks1.shape[0]
        n2 = peaks2.shape[0]
        if n1 == 0 or n2 == 0:
            return 0.0, 0.0

        min_mz = max(peaks2[0, 0], mass_begin)
        max_mz = min(peaks2[n2 - 1, 0], mass_end)
        focused_mz = min_mz
        remaind_m = 0
        remaind_l = 0
        max_lib_intensity = 0.0
        for i in range(n2):
            if peaks2[i, 1] > max_lib_intensity:
                max_lib_intensity = peaks2[i, 1]

        counter = 0
        lib_counter = 0

        while focused_mz <= max_mz:
            sum_l = 0.0
            i = remaind_l
            while i < n2:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_l += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            if sum_l >= 0.01 * max_lib_intensity:
                lib_counter += 1

            sum_m = 0.0
            i = remaind_m
            while i < n1:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            if sum_m > 0 and sum_l >= 0.01 * max_lib_intensity:
                counter += 1

            if focused_mz + bin_width > peaks2[n2 - 1, 0]:
                break
            if remaind_l < n2:
                focused_mz = peaks2[remaind_l, 0]
            else:
                break

        if lib_counter == 0:
            return 0.0, 0.0
        return counter / lib_counter, float(counter)

    @nb.njit
    def _enhanced_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end, penalty):
        n1 = peaks1.shape[0]
        n2 = peaks2.shape[0]
        if n1 == 0 or n2 == 0:
            return 0.0

        min_mz = max(mass_begin, min(peaks1[0, 0], peaks2[0, 0]))
        max_mz = min(mass_end, max(peaks1[n1 - 1, 0], peaks2[n2 - 1, 0]))
        focused_mz = min_mz
        remaind_m = 0
        remaind_l = 0

        max_size = n1 + n2
        mz_arr = np.empty(max_size, dtype=np.float64)
        sum_m_arr = np.empty(max_size, dtype=np.float64)
        sum_r_arr = np.empty(max_size, dtype=np.float64)
        size = 0
        base_m = 0.0
        base_r = 0.0

        while focused_mz <= max_mz:
            sum_m = 0.0
            i = remaind_m
            while i < n1:
                mz = peaks1[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_m += peaks1[i, 1]
                    i += 1
                else:
                    break
            remaind_m = i

            sum_r = 0.0
            i = remaind_l
            while i < n2:
                mz = peaks2[i, 0]
                if mz < focused_mz - bin_width:
                    i += 1
                    continue
                if focused_mz - bin_width <= mz < focused_mz + bin_width:
                    sum_r += peaks2[i, 1]
                    i += 1
                else:
                    break
            remaind_l = i

            mz_arr[size] = focused_mz
            sum_m_arr[size] = sum_m
            sum_r_arr[size] = sum_r
            size += 1
            if sum_m > base_m:
                base_m = sum_m
            if sum_r > base_r:
                base_r = sum_r

            if focused_mz + bin_width > max(peaks1[n1 - 1, 0], peaks2[n2 - 1, 0]):
                break
            next_m = peaks1[remaind_m, 0] if remaind_m < n1 else np.inf
            next_l = peaks2[remaind_l, 0] if remaind_l < n2 else np.inf
            if focused_mz + bin_width > next_l and focused_mz + bin_width <= next_m:
                focused_mz = next_m
            elif focused_mz + bin_width <= next_l and focused_mz + bin_width > next_m:
                focused_mz = next_l
            else:
                focused_mz = min(next_m, next_l)

        if base_m == 0 or base_r == 0:
            return 0.0

        l_counter = 0
        for i in range(size):
            r_norm = sum_r_arr[i] / base_r
            if r_norm > 0.1:
                l_counter += 1

        peak_count_penalty = 1.0
        if l_counter == 1:
            peak_count_penalty = 0.75
        elif l_counter == 2:
            peak_count_penalty = 0.88
        elif l_counter == 3:
            peak_count_penalty = 0.94
        elif l_counter == 4:
            peak_count_penalty = 0.97

        scalar_m = 0.0
        scalar_r = 0.0
        covariance = 0.0
        cutoff = 0.01
        for i in range(size):
            r_int = sum_r_arr[i] / base_r
            if r_int < cutoff:
                continue
            m_int = sum_m_arr[i] / base_m
            mz = mz_arr[i]
            scalar_r += r_int * mz
            if m_int == 0.0:
                scalar_m += m_int * (1.0 - penalty) * mz
            else:
                scalar_m += m_int * mz
            covariance += np.sqrt(m_int * r_int) * mz

        if scalar_m == 0 or scalar_r == 0:
            return 0.0
        return (covariance ** 2) / scalar_m / scalar_r * peak_count_penalty

    # ------------------------------------------------------------------
    #  Matrix-level wrappers (parallel over queries)
    # ------------------------------------------------------------------
    @nb.njit(parallel=True)
    def _simple_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                   bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in nb.prange(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _simple_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end)
        return result

    @nb.njit(parallel=True)
    def _weighted_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                     bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in nb.prange(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _weighted_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end)
        return result

    @nb.njit(parallel=True)
    def _reverse_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                    bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in nb.prange(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _reverse_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end)
        return result

    @nb.njit(parallel=True)
    def _matched_peaks_scores_matrix(peaks, offsets, other_peaks, other_offsets,
                                     bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        ratio = np.empty((m, k), dtype=np.float64)
        count = np.empty((m, k), dtype=np.float64)
        for i in nb.prange(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                r, c = _matched_peaks_scores_pair(p1, p2, bin_width, mass_begin, mass_end)
                ratio[i, j] = r
                count[i, j] = c
        return ratio, count

    @nb.njit(parallel=True)
    def _enhanced_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                     bin_width, mass_begin, mass_end, penalty):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in nb.prange(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _enhanced_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end, penalty)
        return result

else:
    # Pure-Python fallback: re-use the MSdata instance methods.
    def _simple_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        spec1 = MSdata(peaks1, to_normalized=False)
        spec2 = MSdata(peaks2, to_normalized=False)
        return spec1.simple_dot_product(spec2, bin_width, mass_begin, mass_end)

    def _weighted_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        spec1 = MSdata(peaks1, to_normalized=False)
        spec2 = MSdata(peaks2, to_normalized=False)
        return spec1.weighted_dot_product(spec2, bin_width, mass_begin, mass_end)

    def _reverse_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        spec1 = MSdata(peaks1, to_normalized=False)
        spec2 = MSdata(peaks2, to_normalized=False)
        return spec1.reverse_dot_product(spec2, bin_width, mass_begin, mass_end)

    def _matched_peaks_scores_pair(peaks1, peaks2, bin_width, mass_begin, mass_end):
        spec1 = MSdata(peaks1, to_normalized=False)
        spec2 = MSdata(peaks2, to_normalized=False)
        scores = spec1.matched_peaks_scores(spec2, bin_width, mass_begin, mass_end)
        return scores[0], float(scores[1])

    def _enhanced_dot_product_pair(peaks1, peaks2, bin_width, mass_begin, mass_end, penalty):
        spec1 = MSdata(peaks1, to_normalized=False)
        spec2 = MSdata(peaks2, to_normalized=False)
        return spec1.enhanced_dot_product(spec2, bin_width, mass_begin, mass_end, penalty)

    def _simple_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                   bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in range(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _simple_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end)
        return result

    def _weighted_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                     bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in range(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _weighted_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end)
        return result

    def _reverse_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                    bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in range(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _reverse_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end)
        return result

    def _matched_peaks_scores_matrix(peaks, offsets, other_peaks, other_offsets,
                                     bin_width, mass_begin, mass_end):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        ratio = np.empty((m, k), dtype=np.float64)
        count = np.empty((m, k), dtype=np.float64)
        for i in range(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                r, c = _matched_peaks_scores_pair(p1, p2, bin_width, mass_begin, mass_end)
                ratio[i, j] = r
                count[i, j] = c
        return ratio, count

    def _enhanced_dot_product_matrix(peaks, offsets, other_peaks, other_offsets,
                                     bin_width, mass_begin, mass_end, penalty):
        m = offsets.shape[0] - 1
        k = other_offsets.shape[0] - 1
        result = np.empty((m, k), dtype=np.float64)
        for i in range(m):
            p1 = peaks[offsets[i]:offsets[i + 1]]
            for j in range(k):
                p2 = other_peaks[other_offsets[j]:other_offsets[j + 1]]
                result[i, j] = _enhanced_dot_product_pair(p1, p2, bin_width, mass_begin, mass_end, penalty)
        return result




# ---------------------------------------------------------------------------
#  CuPy / GPU availability probe
# ---------------------------------------------------------------------------
try:
    import cupy as _cp
    from cupyx.scipy import sparse as _cupy_sparse
    _HAS_CUPY = True
except Exception:
    _cp = None
    _cupy_sparse = None
    _HAS_CUPY = False


# ---------------------------------------------------------------------------
#  Grid-based (fixed-bin) similarity helpers
#  These provide an approximate but fully-vectorised path for large-scale
#  M × K matrix computation.  When CuPy is available the sparse matrices
#  are kept on the GPU.
# ---------------------------------------------------------------------------

def _grid_prepare(peaks, offsets, bin_width, mass_begin, mass_end):
    """Convert a CSR-like batch to a sparse grid matrix (CPU)."""
    from scipy.sparse import csr_matrix
    n_spectra = offsets.shape[0] - 1
    n_bins = int(np.ceil((mass_end - mass_begin) / bin_width))
    rows = []
    cols = []
    data = []
    for s in range(n_spectra):
        spec = peaks[offsets[s]:offsets[s + 1]]
        for i in range(spec.shape[0]):
            mz = spec[i, 0]
            if mz < mass_begin or mz >= mass_end:
                continue
            bin_idx = int((mz - mass_begin) / bin_width)
            rows.append(s)
            cols.append(bin_idx)
            data.append(spec[i, 1])
    grid = csr_matrix((data, (rows, cols)), shape=(n_spectra, n_bins), dtype=np.float64)
    bin_centers = np.arange(n_bins) * bin_width + mass_begin + bin_width / 2.0
    return grid, bin_centers


def _grid_to_device(grid, bin_centers):
    """Move a CPU sparse grid to the GPU when CuPy is available."""
    if _HAS_CUPY:
        return _cupy_sparse.csr_matrix(grid), _cp.asarray(bin_centers)
    return grid, bin_centers


def _grid_simple_dot_product(grid1, grid2):
    """Approximate simple dot-product via fixed-bin sparse grids."""
    xp = _cp if _HAS_CUPY and hasattr(grid1, 'device') else np
    scipy_sparse = __import__('scipy.sparse', fromlist=[''])

    # base-peak normalisation to 999 (matches exact MSDIAL logic)
    max1 = grid1.max(axis=1).toarray().ravel()
    max2 = grid2.max(axis=1).toarray().ravel()
    max1[max1 == 0] = 1.0
    max2[max2 == 0] = 1.0

    if xp is _cp:
        inv_max1 = _cp.asarray(1.0 / max1)
        inv_max2 = _cp.asarray(1.0 / max2)
        g1 = grid1.multiply(inv_max1[:, None]).multiply(999.0).tocsr()
        g2 = grid2.multiply(inv_max2[:, None]).multiply(999.0).tocsr()
        g1_sqrt = _cupy_sparse.csr_matrix(
            (xp.sqrt(g1.data), g1.indices, g1.indptr), shape=g1.shape
        )
        g2_sqrt = _cupy_sparse.csr_matrix(
            (xp.sqrt(g2.data), g2.indices, g2.indptr), shape=g2.shape
        )
        cov = g1_sqrt.dot(g2_sqrt.T).toarray()
        scalar_m = xp.asarray(g1.sum(axis=1)).ravel()[:, None]
        scalar_r = xp.asarray(g2.sum(axis=1)).ravel()[None, :]
    else:
        g1 = grid1.multiply((1.0 / max1)[:, None]).multiply(999.0).tocsr()
        g2 = grid2.multiply((1.0 / max2)[:, None]).multiply(999.0).tocsr()
        g1_sqrt = scipy_sparse.csr_matrix(
            (np.sqrt(g1.data), g1.indices, g1.indptr), shape=g1.shape
        )
        g2_sqrt = scipy_sparse.csr_matrix(
            (np.sqrt(g2.data), g2.indices, g2.indptr), shape=g2.shape
        )
        cov = g1_sqrt.dot(g2_sqrt.T).toarray()
        scalar_m = np.asarray(g1.sum(axis=1)).ravel()[:, None]
        scalar_r = np.asarray(g2.sum(axis=1)).ravel()[None, :]

    with np.errstate(divide='ignore', invalid='ignore'):
        result = (cov ** 2) / (scalar_m * scalar_r)
    if xp is _cp:
        result = _cp.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
        return _cp.asnumpy(result)
    return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)


def _grid_weighted_or_reverse_dot_product(grid1, grid2, bin_centers, reverse=False, enhanced=False, penalty=0.5):
    """Approximate weighted / reverse / enhanced dot-product via fixed-bin sparse grids."""
    xp = _cp if _HAS_CUPY and hasattr(grid1, 'device') else np
    scipy_sparse = __import__('scipy.sparse', fromlist=[''])

    # normalise
    max1 = grid1.max(axis=1).toarray().ravel()
    max2 = grid2.max(axis=1).toarray().ravel()
    max1[max1 == 0] = 1.0
    max2[max2 == 0] = 1.0

    if xp is _cp:
        inv_max1 = _cp.asarray(1.0 / max1)
        inv_max2 = _cp.asarray(1.0 / max2)
        g1 = grid1.multiply(inv_max1[:, None]).tocsr()
        g2 = grid2.multiply(inv_max2[:, None]).tocsr()
    else:
        g1 = grid1.multiply((1.0 / max1)[:, None]).tocsr()
        g2 = grid2.multiply((1.0 / max2)[:, None]).tocsr()

    # peak-count penalty from reference (g2)
    l_counter = np.asarray((g2 > 0.1).sum(axis=1)).ravel().astype(np.int64)
    penalty_vec = np.ones(g2.shape[0], dtype=np.float64)
    for i, c in enumerate(l_counter):
        if c == 1:
            penalty_vec[i] = 0.75
        elif c == 2:
            penalty_vec[i] = 0.88
        elif c == 3:
            penalty_vec[i] = 0.94
        elif c == 4:
            penalty_vec[i] = 0.97
    if xp is _cp:
        penalty_vec = _cp.asarray(penalty_vec)

    # apply cutoff
    cutoff = 0.01
    if reverse or enhanced:
        # cutoff on reference (g2)
        g2_cut = g2.copy()
        if xp is _cp:
            g2_cut.data[g2_cut.data < cutoff] = 0.0
        else:
            g2_cut.data[g2_cut.data < cutoff] = 0.0
        g2_cut.eliminate_zeros()
        g1_cut = g1
    else:
        # cutoff on measured (g1)
        g1_cut = g1.copy()
        if xp is _cp:
            g1_cut.data[g1_cut.data < cutoff] = 0.0
        else:
            g1_cut.data[g1_cut.data < cutoff] = 0.0
        g1_cut.eliminate_zeros()
        g2_cut = g2

    # weighted vectors: sqrt(intensity) * sqrt(mz)
    sqrt_mz = xp.sqrt(np.asarray(bin_centers) if xp is not _cp else bin_centers)
    if xp is _cp:
        d_sqrt_mz = _cupy_sparse.diags(sqrt_mz)
    else:
        d_sqrt_mz = scipy_sparse.diags(sqrt_mz)

    # sqrt of non-zero data
    if xp is _cp:
        g1_sqrt = _cupy_sparse.csr_matrix(
            (xp.sqrt(g1_cut.data), g1_cut.indices, g1_cut.indptr), shape=g1_cut.shape
        )
        g2_sqrt = _cupy_sparse.csr_matrix(
            (xp.sqrt(g2_cut.data), g2_cut.indices, g2_cut.indptr), shape=g2_cut.shape
        )
    else:
        g1_sqrt = scipy_sparse.csr_matrix(
            (np.sqrt(g1_cut.data), g1_cut.indices, g1_cut.indptr), shape=g1_cut.shape
        )
        g2_sqrt = scipy_sparse.csr_matrix(
            (np.sqrt(g2_cut.data), g2_cut.indices, g2_cut.indptr), shape=g2_cut.shape
        )

    g1_w = g1_sqrt.dot(d_sqrt_mz)
    g2_w = g2_sqrt.dot(d_sqrt_mz)

    # covariance
    cov = g1_w.dot(g2_w.T).toarray()

    # scalarM & scalarR via original intensities (not sqrt) weighted by mz
    if xp is _cp:
        d_mz = _cupy_sparse.diags(np.asarray(bin_centers) if xp is not _cp else bin_centers)
    else:
        d_mz = scipy_sparse.diags(np.asarray(bin_centers) if xp is not _cp else bin_centers)

    if reverse or enhanced:
        # cutoff is on reference (g2).  scalar_r already matches exact because
        # g2_cut contains only bins above cutoff.  scalar_m must be restricted
        # to the same columns so that unmatched measured bins are excluded.
        g2_bool = (g2_cut > 0).astype(np.float64)
        if xp is _cp:
            g2_bool = _cupy_sparse.csr_matrix(g2_bool)
        W = g1_cut.dot(d_mz)
        scalar_m = W.dot(g2_bool.T).toarray()
        scalar_r = np.asarray(g2_cut.dot(d_mz).sum(axis=1)).ravel()[None, :]
    else:
        # cutoff is on measured (g1).  scalar_m already matches exact because
        # g1_cut contains only bins above cutoff.  scalar_r must be restricted
        # to the same columns so that reference peaks opposite a filtered
        # measured bin are excluded.
        g1_bool = (g1_cut > 0).astype(np.float64)
        if xp is _cp:
            g1_bool = _cupy_sparse.csr_matrix(g1_bool)
        W = g2.dot(d_mz)
        scalar_r = g1_bool.dot(W.T).toarray()
        scalar_m = np.asarray(g1_cut.dot(d_mz).sum(axis=1)).ravel()[:, None]

    with np.errstate(divide='ignore', invalid='ignore'):
        result = (cov ** 2) / (scalar_m * scalar_r)
    result = result * penalty_vec[None, :]

    if xp is _cp:
        result = _cp.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
        return _cp.asnumpy(result)
    return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)


def _grid_matched_peaks_scores(grid1, grid2):
    """Approximate matched-peaks scores via fixed-bin sparse grids."""
    xp = _cp if _HAS_CUPY and hasattr(grid1, 'device') else np

    # max intensity per reference spectrum
    max2 = grid2.max(axis=1).toarray().ravel()
    max2[max2 == 0] = 1.0

    if xp is _cp:
        threshold = _cp.asarray(max2 * 0.01)
        g2_thresh = grid2.multiply(threshold[:, None] ** -1) >= 1.0
        # bool sparse -> int sparse
        g2_int = g2_thresh.astype(np.float64)
        g1_bool = (grid1 > 0).astype(np.float64)
        lib_counter = g2_int.sum(axis=1).ravel()
        matched = g1_bool.dot(g2_int.T).toarray()
    else:
        g2_norm = grid2.multiply((1.0 / max2)[:, None])
        g2_thresh = (g2_norm >= 0.01).astype(np.float64)
        g1_bool = (grid1 > 0).astype(np.float64)
        lib_counter = np.array(g2_thresh.sum(axis=1)).ravel()
        matched = g1_bool.dot(g2_thresh.T).toarray()

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = matched / lib_counter[None, :]
    if xp is _cp:
        ratio = _cp.nan_to_num(ratio, nan=0.0, posinf=0.0, neginf=0.0)
        return _cp.asnumpy(ratio), _cp.asnumpy(matched)
    ratio = np.nan_to_num(ratio, nan=0.0, posinf=0.0, neginf=0.0)
    return ratio, matched



class MSBatch:
    """
    Compact CSR-like container for a collection of MS spectra.

    ``MSBatch`` stores spectra in a flattened ``peaks`` array together with
    ``offsets`` (cumulated lengths).  This eliminates zero-padding overhead
    and works efficiently with the Numba-accelerated matrix kernels.

    Two compute backends are supported:

    * ``device='cpu'`` (default) – exact MS-DIAL dynamic-binning algorithm
      accelerated by Numba (or a pure-Python fallback).
    * ``device='cuda'`` – fixed-grid projection onto a regular m/z lattice.
      The grids are built as sparse matrices and all algebra is performed on
      the GPU via CuPy when available, otherwise it falls back to CPU sparse
      algebra (still much faster than pairwise Python loops for large batches).

    Parameters
    ----------
    peaks : np.ndarray, shape (total_peaks, 2)
        Flattened (mz, intensity) array.
    offsets : np.ndarray, shape (n_spectra + 1,)
        Cumulative offsets; spectrum *i* occupies ``peaks[offsets[i]:offsets[i+1]]``.
    device : {'cpu', 'cuda'}, default 'cpu'
        Compute backend.  ``'cuda'`` requires CuPy; if it is missing the
        instance silently falls back to ``'cpu'`` and issues a warning.

    Examples
    --------
    >>> from mzpy.ms import MSdata, MSBatch
    >>> spec1 = MSdata([[100.0, 50], [200.0, 100]], to_normalized=False)
    >>> spec2 = MSdata([[150.0, 60]], to_normalized=False)
    >>> batch = MSBatch.from_list([spec1, spec2])
    >>> len(batch)
    2
    >>> batch.get_spectrum(0)
    array([[100.,  50.],
           [200., 100.]])
    """

    def __init__(self, peaks: np.ndarray, offsets: np.ndarray, device: str = 'cpu'):
        self.peaks = np.asarray(peaks, dtype=np.float64)
        self.offsets = np.asarray(offsets, dtype=np.int64)
        if self.offsets[0] != 0:
            raise ValueError("offsets must start with 0")
        if self.offsets[-1] != self.peaks.shape[0]:
            raise ValueError("offsets[-1] must equal total number of peaks")
        self.n_spectra = self.offsets.shape[0] - 1

        dev = device.lower()
        if dev == 'cuda' and not _HAS_CUPY:
            import warnings
            warnings.warn(
                "CuPy is not installed or GPU is not available; "
                "grid path will use CPU sparse algebra for MSBatch."
            )
            # Keep dev == 'cuda' so _use_grid_path still returns True;
            # the grid helpers already fall back to CPU sparse algebra.
        self.device = dev

        # cached grid data (fixed-bin projection)
        self._grid = None
        self._bin_centers = None
        self._grid_params = None

    # ------------------------------------------------------------------
    #  Construction helpers
    # ------------------------------------------------------------------
    @classmethod
    def from_list(cls, spectra: Sequence[Union["MSdata", List[List[float]], np.ndarray]], device: str = 'cpu'):
        """
        Build an ``MSBatch`` from a list of spectra.

        Each element can be an ``MSdata`` instance, a nested Python list
        ``[[mz, int], ...]``, or a ``(N, 2)`` ndarray.  Peaks with
        non-positive m/z are discarded and each spectrum is sorted by m/z.
        """
        arrays = []
        for s in spectra:
            arr = np.asarray(s)
            if arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError(
                    f"Each spectrum must have shape (N, 2); got {arr.shape}"
                )
            # keep only rows with positive m/z
            arr = arr[arr[:, 0] > 0]
            if arr.shape[0] > 1:
                arr = arr[arr[:, 0].argsort()]
            arrays.append(arr)

        total = sum(a.shape[0] for a in arrays)
        peaks = np.empty((total, 2), dtype=np.float64)
        offsets = np.empty(len(arrays) + 1, dtype=np.int64)
        offsets[0] = 0
        pos = 0
        for i, a in enumerate(arrays):
            peaks[pos:pos + a.shape[0]] = a
            pos += a.shape[0]
            offsets[i + 1] = pos
        return cls(peaks, offsets, device=device)

    @classmethod
    def from_stack(cls, stacked: np.ndarray, device: str = 'cpu'):
        """
        Build an ``MSBatch`` from a zero-padded ``(M, N, 2)`` stack array
        (as produced by ``MSdata.stack()``).

        Rows with non-positive m/z are treated as padding and stripped.
        """
        stacked = np.asarray(stacked)
        if stacked.ndim != 3 or stacked.shape[2] != 2:
            raise ValueError(f"stacked must have shape (M, N, 2); got {stacked.shape}")
        m, n, _ = stacked.shape
        lengths = np.count_nonzero(stacked[:, :, 0] > 0, axis=1)
        total = int(lengths.sum())
        peaks = np.empty((total, 2), dtype=np.float64)
        offsets = np.empty(m + 1, dtype=np.int64)
        offsets[0] = 0
        pos = 0
        for i in range(m):
            l = int(lengths[i])
            arr = stacked[i, :l]
            if l > 1:
                arr = arr[arr[:, 0].argsort()]
            peaks[pos:pos + l] = arr
            pos += l
            offsets[i + 1] = pos
        return cls(peaks, offsets, device=device)

    def __len__(self):
        return self.n_spectra

    def __repr__(self):
        return f"MSBatch(n_spectra={self.n_spectra}, total_peaks={self.peaks.shape[0]}, device='{self.device}')"

    def get_spectrum(self, idx: int) -> np.ndarray:
        """Return the *idx*-th spectrum as a ``(N, 2)`` ndarray."""
        if idx < 0 or idx >= self.n_spectra:
            raise IndexError(f"index {idx} out of range for {self.n_spectra} spectra")
        return self.peaks[self.offsets[idx]:self.offsets[idx + 1]]

    def to_stack(self, n: int = None) -> np.ndarray:
        """
        Convert back to a zero-padded ``(M, N, 2)`` array.

        Parameters
        ----------
        n : int, optional
            Fixed padding length.  If None, uses the maximum spectrum length.
        """
        if n is None:
            n = max(self.offsets[i + 1] - self.offsets[i] for i in range(self.n_spectra))
        result = np.zeros((self.n_spectra, n, 2), dtype=np.float64)
        for i in range(self.n_spectra):
            spec = self.get_spectrum(i)
            k = min(spec.shape[0], n)
            result[i, :k] = spec[:k]
        return result

    # ------------------------------------------------------------------
    #  Grid preparation (fixed-bin projection)
    # ------------------------------------------------------------------
    def prepare_grid(self, bin: float = 0.05, mass_begin: float = 0.0, mass_end: float = 2000.0):
        """
        Project all spectra onto a regular m/z grid.

        Peaks falling inside the same bin are summed.  The resulting sparse
        matrix is cached inside the batch instance and reused by subsequent
        matrix methods when ``device='cuda'`` (or when both batches use the
        grid path).

        Returns
        -------
        grid : scipy.sparse.csr_matrix or cupyx.scipy.sparse.csr_matrix
        bin_centers : np.ndarray or cupy.ndarray
        """
        grid_cpu, bin_centers_cpu = _grid_prepare(
            self.peaks, self.offsets, bin, mass_begin, mass_end
        )
        if self.device == 'cuda' and _HAS_CUPY:
            self._grid, self._bin_centers = _grid_to_device(grid_cpu, bin_centers_cpu)
        else:
            self._grid, self._bin_centers = grid_cpu, bin_centers_cpu
        self._grid_params = (bin, mass_begin, mass_end)
        return self._grid, self._bin_centers

    def _ensure_grid(self, other, bin, mass_begin, mass_end):
        """Internal helper: build grids for both sides when needed."""
        if self._grid is None or self._grid_params != (bin, mass_begin, mass_end):
            self.prepare_grid(bin, mass_begin, mass_end)
        if other._grid is None or other._grid_params != (bin, mass_begin, mass_end):
            other.prepare_grid(bin, mass_begin, mass_end)

    def _use_grid_path(self, other):
        """True if either side requests the grid (cuda) path."""
        return self.device == 'cuda' or other.device == 'cuda'

    # ------------------------------------------------------------------
    #  Batch similarity matrices
    # ------------------------------------------------------------------
    def simple_dot_product_matrix(self, other: "MSBatch", bin: float = 0.05,
                                  mass_begin: float = 0.0, mass_end: float = 2000.0) -> np.ndarray:
        """
        Compute the simple dot-product similarity matrix.

        Returns
        -------
        np.ndarray, shape (len(self), len(other))
        """
        if self._use_grid_path(other):
            self._ensure_grid(other, bin, mass_begin, mass_end)
            return _grid_simple_dot_product(self._grid, other._grid)
        return _simple_dot_product_matrix(
            self.peaks, self.offsets, other.peaks, other.offsets,
            bin, mass_begin, mass_end
        )

    def weighted_dot_product_matrix(self, other: "MSBatch", bin: float = 0.05,
                                    mass_begin: float = 0.0, mass_end: float = 2000.0) -> np.ndarray:
        """
        Compute the weighted dot-product similarity matrix.

        Returns
        -------
        np.ndarray, shape (len(self), len(other))
        """
        if self._use_grid_path(other):
            self._ensure_grid(other, bin, mass_begin, mass_end)
            return _grid_weighted_or_reverse_dot_product(
                self._grid, other._grid, self._bin_centers, reverse=False, enhanced=False
            )
        return _weighted_dot_product_matrix(
            self.peaks, self.offsets, other.peaks, other.offsets,
            bin, mass_begin, mass_end
        )

    def reverse_dot_product_matrix(self, other: "MSBatch", bin: float = 0.05,
                                   mass_begin: float = 0.0, mass_end: float = 2000.0) -> np.ndarray:
        """
        Compute the reverse dot-product similarity matrix.

        Returns
        -------
        np.ndarray, shape (len(self), len(other))
        """
        if self._use_grid_path(other):
            self._ensure_grid(other, bin, mass_begin, mass_end)
            return _grid_weighted_or_reverse_dot_product(
                self._grid, other._grid, self._bin_centers, reverse=True, enhanced=False
            )
        return _reverse_dot_product_matrix(
            self.peaks, self.offsets, other.peaks, other.offsets,
            bin, mass_begin, mass_end
        )

    def matched_peaks_scores_matrix(self, other: "MSBatch", bin: float = 0.05,
                                    mass_begin: float = 0.0, mass_end: float = 2000.0):
        """
        Compute the matched-peaks scores matrix.

        Returns
        -------
        ratio : np.ndarray, shape (len(self), len(other))
        count : np.ndarray, shape (len(self), len(other))
        """
        if self._use_grid_path(other):
            self._ensure_grid(other, bin, mass_begin, mass_end)
            ratio, count = _grid_matched_peaks_scores(self._grid, other._grid)
            return ratio, count
        ratio, count = _matched_peaks_scores_matrix(
            self.peaks, self.offsets, other.peaks, other.offsets,
            bin, mass_begin, mass_end
        )
        return ratio, count

    def enhanced_dot_product_matrix(self, other: "MSBatch", bin: float = 0.05,
                                    mass_begin: float = 0.0, mass_end: float = 2000.0,
                                    penalty: float = 0.5) -> np.ndarray:
        """
        Compute the enhanced dot-product similarity matrix.

        Returns
        -------
        np.ndarray, shape (len(self), len(other))
        """
        if self._use_grid_path(other):
            self._ensure_grid(other, bin, mass_begin, mass_end)
            return _grid_weighted_or_reverse_dot_product(
                self._grid, other._grid, self._bin_centers,
                reverse=False, enhanced=True, penalty=penalty
            )
        return _enhanced_dot_product_matrix(
            self.peaks, self.offsets, other.peaks, other.offsets,
            bin, mass_begin, mass_end, penalty
        )


# ---------------------------------------------------------------------------
#  Example / quick-smoke test when run as a script
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import warnings
    import time

    print("=" * 60)
    print("mzpy.ms – quick smoke test")
    print("=" * 60)

    # --------------------------------------------------------------
    #  1. MSdata basics
    # --------------------------------------------------------------
    print("\n1. MSdata pairwise similarity")
    spec_q = MSdata([[100.0, 50], [200.0, 100], [300.0, 30]], to_normalized=False)
    spec_r = MSdata([[100.0, 60], [201.0, 90]], to_normalized=False)

    print(f"   query : {spec_q}")
    print(f"   ref   : {spec_r}")
    print(f"   simple_dot_product      = {spec_q.simple_dot_product(spec_r):.4f}")
    print(f"   weighted_dot_product    = {spec_q.weighted_dot_product(spec_r):.4f}")
    print(f"   reverse_dot_product     = {spec_q.reverse_dot_product(spec_r):.4f}")
    print(f"   enhanced_dot_product    = {spec_q.enhanced_dot_product(spec_r):.4f}")
    ratio, count = spec_q.matched_peaks_scores(spec_r)
    print(f"   matched_peaks_scores    = ratio {ratio:.4f}, count {count}")
    print(f"   gaussian_similarity     = {spec_q.gaussian_similarity(spec_r):.4f}")
    print(f"   spectral_entropy_sim    = {spec_q.spectral_entropy_similarity(spec_r):.4f}")

    # --------------------------------------------------------------
    #  2. MSBatch – exact CPU path (default device='cpu')
    # --------------------------------------------------------------
    print("\n2. MSBatch exact CPU path")
    a = MSdata([[100.0, 1.0], [200.0, 2.0]], to_normalized=False)
    b = MSdata([[100.0, 1.0], [200.0, 3.0]], to_normalized=False)
    c = MSdata([[100.0, 2.0], [300.0, 1.0]], to_normalized=False)

    batch_q = MSBatch.from_list([a, b])
    batch_r = MSBatch.from_list([a, c])

    sdp_cpu = batch_q.simple_dot_product_matrix(batch_r)
    wdp_cpu = batch_q.weighted_dot_product_matrix(batch_r)
    rdp_cpu = batch_q.reverse_dot_product_matrix(batch_r)
    edp_cpu = batch_q.enhanced_dot_product_matrix(batch_r)
    mp_cpu = batch_q.matched_peaks_scores_matrix(batch_r)

    print(f"   batch_q: {batch_q}")
    print(f"   batch_r: {batch_r}")
    print(f"   simple_dot_product_matrix:\n{sdp_cpu}")
    print(f"   weighted_dot_product_matrix:\n{wdp_cpu}")
    print(f"   reverse_dot_product_matrix:\n{rdp_cpu}")
    print(f"   matched_peaks ratio:\n{mp_cpu[0]}")

    # --------------------------------------------------------------
    #  3. MSBatch – grid path (device='cuda' falls back to CPU sparse)
    # --------------------------------------------------------------
    print("\n3. MSBatch grid path (device='cuda' → CPU sparse fallback)")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        batch_qg = MSBatch.from_list([a, b], device="cuda")
        batch_rg = MSBatch.from_list([a, c], device="cuda")
        if w:
            print(f"   (warning) {w[-1].message}")

    sdp_grid = batch_qg.simple_dot_product_matrix(batch_rg)
    wdp_grid = batch_qg.weighted_dot_product_matrix(batch_rg)
    rdp_grid = batch_qg.reverse_dot_product_matrix(batch_rg)
    edp_grid = batch_qg.enhanced_dot_product_matrix(batch_rg)
    mp_grid = batch_qg.matched_peaks_scores_matrix(batch_rg)

    print(f"   simple max_diff  = {np.max(np.abs(sdp_cpu - sdp_grid)):.3e}")
    print(f"   weighted max_diff= {np.max(np.abs(wdp_cpu - wdp_grid)):.3e}")
    print(f"   reverse max_diff = {np.max(np.abs(rdp_cpu - rdp_grid)):.3e}")
    print(f"   matched max_diff = {np.max(np.abs(mp_cpu[0] - mp_grid[0])):.3e}")

    # --------------------------------------------------------------
    #  4. Performance smoke test (random spectra)
    # --------------------------------------------------------------
    print("\n4. Performance smoke test (100 × 100 spectra, ~75 peaks each)")
    rng = np.random.default_rng(42)
    specs_q = []
    specs_r = []
    for _ in range(100):
        n = rng.integers(50, 100)
        mz = np.sort(rng.uniform(50, 1000, n))
        inten = rng.exponential(100, n)
        specs_q.append(MSdata(np.column_stack([mz, inten]), to_normalized=False))
        n = rng.integers(50, 100)
        mz = np.sort(rng.uniform(50, 1000, n))
        inten = rng.exponential(100, n)
        specs_r.append(MSdata(np.column_stack([mz, inten]), to_normalized=False))

    bq = MSBatch.from_list(specs_q)
    br = MSBatch.from_list(specs_r)
    bqg = MSBatch.from_list(specs_q, device="cuda")
    brg = MSBatch.from_list(specs_r, device="cuda")

    t0 = time.time()
    _ = bq.simple_dot_product_matrix(br)
    t_cpu = time.time() - t0

    t0 = time.time()
    _ = bqg.simple_dot_product_matrix(brg)
    t_grid = time.time() - t0

    print(f"   CPU exact : {t_cpu:.3f} s")
    print(f"   Grid      : {t_grid:.3f} s")
    print(f"   speed-up  : {t_cpu / t_grid:.1f}×")

    print("\n" + "=" * 60)
    print("All smoke tests finished.")
    print("=" * 60)
