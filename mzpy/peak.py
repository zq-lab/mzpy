# -*- coding: utf-8 -*-
'''
更新说明：
    剔除了所有的绘图代码，轻量化
    绘图应交给用户调用自己喜欢的绘图包处理，或者使用mzpy独立的绘图模块
    绝大多数情况下，seaborn搭配数据框已足够好用
    同时有利于降低耦合度
'''
'''
MSMS database or data sheet processor

MSP:
    In every filed, it cannot contain newline character, such as "\r" or "\n".
        Or the msp text structure will be interrupted after being output.
        Therefore, in the PeakFrame.to_msp method, "\r" and "\n" were checked and deleted firstly.
    FORMULA can be ''. But it can not be "nan" which can not be accpted by MS-Dial.
    MS-Dial does not accept single autom or ion, such as Na, N, S. 
        Thus items also be checked atom bumber before being exported in to_msp function.
'''

import ast
import numpy as np
import pandas as pd
import re
import warnings

from . import ms, mz
from . import similarity
from .rest import np_classfy_df
from .stat import enrich_df


__FIELDS__ = (
# Standard Fields for LC-MS Ion Peak Information
#   Num Peaks was NOT included because this field always needs 
#       to be automatically updated when used. 
#   It is a field that needs to refresh automatically at all times.
              'NAME',           # compound name
              'PRECURSORMZ',    # precursoe (MS1 ion) mz
              'PRECURSORTYPE',  # precursor (MS1 ion) type
              'IONMODE',        # ion mode, pos or neg
              'RETENTIONTIME',  # retention time in minutes
              'FORMULA',        # compound formula
              'ONTOLOGY',       # Commonly used classfication of compound structures
              'SMILES',         # compound structure
              'INCHIKEY',       # compound inchikey
              'INSTRUMENTTYPE', # instrument type, QQQ, Q-TOF, Q-Orbitrap
              'INSTRUMENT',     # instrument model              
              'COLLISIONENERGY',# voltage (V) for ion fragment
              'CCS',            # Collision Cross Section, not necessary
              'COMMENT',        # other information, not necessary
              'MSMS'            # fragment ion list: m/z and intensity
)


class PeakFrame(pd.DataFrame):
    '''
    pandas-like
    MSMS handling
    Two fields are fixed: "Num Peaks" and "MSMS"
    '''

    @property
    def _constructor(self):
        '''
        This is to ensure that existing subtypes are preserved even after 
            using pandas functions like merge and concat.
        '''
        return self.__class__
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    @classmethod
    def _parse_msp_txt(cls, txt, sep_ms2= '\t') -> dict:
        data = {}
        data['MSMS'] = []
        lines = txt.strip().split('\n')
        for line in lines:
            if ':' in line:
                d = line.split(':', 1)
                data[d[0]] = d[1].strip() 
            elif line.strip():
                d = line.strip().split(sep_ms2, 1)
                data['MSMS'].append([float(d[0]), float(d[1])])
        return data
    
    @classmethod
    def _parse_cfmid_txt(cls, txt):
        ion = {}
        MSMS = []
        txt = txt.splitlines()
        for line in txt:
            if line.startswith('#In-silico'):
                x = line.split(' ')
                ion['instrumenttype'] = x[1]
                ion['precursortype'] = x[2]
                if x[2][-1] == '-':
                    ion['ionmode'] = 'Negative'
                elif x[2][-1] == '+':
                    ion['ionmode'] = 'Positive'
            elif line.startswith('#PREDICTED'):
                ion['comment'] = line.split(' ',1)[1]
            elif line.startswith('#ID'):
                ion['name'] = line.split('=', 1)[1]
            elif line.startswith('#SMILES'):
                ion['smiles'] = line.split('=', 1)[1]
            elif line.startswith('#InChiKey'):
                ion['inchikey'] = line.split('=', 1)[1]
            elif line.startswith('#Formula'):
                ion['formula'] = line.split('=', 1)[1]
            elif line.startswith('#PMass'):
                ion['precursormz'] = line.split('=', 1)[1]
            elif line.strip() == '': # Prevent the cursor from overflowing in the next check
                ion['collisionenergy'] = '10, 20 40 V'                
                breaknpeaks_on
            elif line[0].isdigit():
                mz, intensity = line.split(' ')[0:2]
                MSMS.append([float(mz), float(intensity)])
        ion['MSMS'] = MSMS
        return ion
    
    def centroid_msms(self):
        def centroid(data):
            return ms.MSdata(data)
        self['MSMS'] = self['MSMS'].apply(lambda x: centroid(x))

#################################################################

    def convert_ms_array(self, ms_on='MSMS', T=False, inplace=False, ms_len=None, npeaks_on='Num Peaks'):
        '''
        Convert MSMS column to numpy float array with optional fixed length
        
        Parameters:
        -----------
        ms_on : str, default='MSMS'
            Column name containing MSMS data
            
        T : bool, default=False
            If True, transpose array to shape (2, N) where:
                - First element: array of m/z values
                - Second element: array of intensity values
            If False, return shape (N, 2) where each row is [mz, intensity]
        
        inplace : bool, default=False
            If True, modify the PeakFrame in-place and return None
            If False, return the converted data without modifying the original
        
        ms_len : int or None, default=None
            Specify fixed length for MSMS arrays:
            - None: keep original length (no padding)
            - 0: use maximum value in npeaks_on column
            - positive int: use this value as fixed length
                If smaller than max peaks count, raise warning and truncate data
            
            When ms_len is set, insufficient peaks are padded with [0, 0] 
            (or [0] for each m/z and intensity if T=True)
        
        npeaks_on : str, default='Num Peaks'
            Column name containing number of peaks, used when ms_len=0
        
        Returns:
        --------
        pandas.DataFrame or None
            - If inplace=False: new DataFrame with MSMS column converted
            - If inplace=True: None (modifies self in-place)
        
        Raises:
        -------
        ValueError
            - If ms_on column not found
            - If npeaks_on column not found when ms_len=0
        
        Warnings:
        ---------
        UserWarning
            - If ms_len is smaller than max peaks count (data will be truncated)
        
        Examples:
        ---------
        >>> pf = PeakFrame(...)
        
        # Variable length arrays
        >>> pf_new = pf.convert_ms_array()  # shape (N, 2) for each row
        >>> pf_new = pf.convert_ms_array(T=True)  # shape (2, N) for each row
        
        # Fixed length arrays (pad to max peaks count)
        >>> pf_new = pf.convert_ms_array(ms_len=0)  # pad to max in 'Num Peaks'
        >>> pf_new = pf.convert_ms_array(ms_len=0, T=True)
        
        # Fixed length arrays with specific length
        >>> pf_new = pf.convert_ms_array(ms_len=50)  # pad to 50 peaks
        >>> pf_new = pf.convert_ms_array(ms_len=50, T=True)
        
        # Modify in-place
        >>> pf.convert_ms_array(ms_len=0, inplace=True)
        '''
        if ms_on not in self.columns:
            raise ValueError(f"Column '{ms_on}' not found in PeakFrame")
        
        # 确定固定长度
        fixed_len = None
        if ms_len is not None:
            if ms_len == 0:
                # 使用"Num Peaks"列的最大值
                if npeaks_on not in self.columns:
                    raise ValueError(f"Column '{npeaks_on}' not found in PeakFrame")
                fixed_len = int(self[npeaks_on].max())
            elif isinstance(ms_len, int) and ms_len > 0:
                # 检查是否小于最大峰数
                if npeaks_on in self.columns:
                    max_peaks = int(self[npeaks_on].max())
                    if ms_len < max_peaks:
                        warnings.warn(
                            f"ms_len={ms_len} is smaller than max peaks count {max_peaks}. "
                            f"Data will be truncated to {ms_len} peaks.",
                            UserWarning,
                            stacklevel=2
                        )
                fixed_len = ms_len
            else:
                raise ValueError(f"ms_len must be None, 0, or a positive integer, got {ms_len}")
        
        # 获取MSMS数据
        msms_data = self[ms_on].values
        
        # 转换为numpy数组
        def convert_spectrum(spectrum):
            """转换单个谱图为numpy数组"""
            arr = np.asarray(spectrum, dtype=np.float64)
            
            # 确保是二维数组 (N, 2)
            if arr.ndim == 1:
                if len(arr) % 2 != 0:
                    raise ValueError(f"Invalid spectrum shape: {arr.shape}")
                arr = arr.reshape(-1, 2)
            elif arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError(f"Spectrum must have shape (N, 2), got {arr.shape}")
            
            # 如果指定了固定长度，进行padding或truncation
            if fixed_len is not None:
                n_peaks = arr.shape[0]
                if n_peaks < fixed_len:
                    # Pad with [0, 0]
                    padding = np.zeros((fixed_len - n_peaks, 2), dtype=np.float64)
                    arr = np.vstack([arr, padding])
                elif n_peaks > fixed_len:
                    # Truncate
                    arr = arr[:fixed_len]
            
            # 根据T参数转置
            if T:
                return arr.T  # shape (2, N)
            else:
                return arr  # shape (N, 2)
        
        # 转换所有数据为列表（避免numpy数组形状不一致的问题）
        converted = [convert_spectrum(x) for x in msms_data]
        
        if inplace:
            # 就地修改
            self[ms_on] = converted
            return None
        else:
            # 返回新数据框
            result = self.copy()
            result[ms_on] = converted
            return result

    def refresh_npeaks(self, npeaks_on='Num Peaks', T=False, ms_on='MSMS'):
        '''
        Refresh (update) the "Num Peaks" column based on MSMS data
        
        Parameters:
        -----------
        npeaks_on : str, default='Num Peaks'
            Column name to store the number of peaks
            
        T : bool, default=False
            Indicates whether MSMS data is transposed:
            - If T=False: MSMS has shape (N, 2), count peaks from first dimension (N)
            - If T=True: MSMS has shape (2, N), count non-zero values in first element (m/z array)
        
        ms_on : str, default='MSMS'
            Column name containing MSMS data
        
        Returns:
        --------
        None
            Modifies npeaks_on column in-place
        
        Examples:
        ---------
        >>> pf = PeakFrame(...)
        
        # For non-transposed data (standard format)
        >>> pf.refresh_npeaks()
        >>> pf.refresh_npeaks(T=False)
        
        # For transposed data
        >>> pf.refresh_npeaks(T=True)
        
        # With custom column names
        >>> pf.refresh_npeaks(npeaks_on='Peak_Count', T=True, ms_on='MS2')
        '''
        if ms_on not in self.columns:
            raise ValueError(f"Column '{ms_on}' not found in PeakFrame")
        
        def count_peaks(spectrum):
            """Count peaks based on data format (transposed or not)"""
            if spectrum is None or (isinstance(spectrum, float) and np.isnan(spectrum)):
                return 0
            
            arr = np.asarray(spectrum)
            
            if arr.size == 0:
                return 0
            
            if T:
                # Transposed format: shape (2, N)
                # Count non-zero values in the first element (m/z array)
                # or count the length if first element is the mz array
                if arr.ndim == 2 and arr.shape[0] == 2:
                    # Count non-zero m/z values
                    n_peaks = np.count_nonzero(arr[0])
                    return int(n_peaks) if n_peaks > 0 else 0
                else:
                    # Fallback: if shape doesn't match, return length
                    return arr.shape[0] if arr.ndim == 1 else arr.shape[-1]
            else:
                # Non-transposed format: shape (N, 2)
                # Count peaks from first dimension
                if arr.ndim == 2:
                    return arr.shape[0]
                elif arr.ndim == 1:
                    # If 1D array, assume it's a flat list of [mz1, inten1, mz2, inten2, ...]
                    return arr.shape[0] // 2 if arr.shape[0] % 2 == 0 else 0
                else:
                    return 0
        
        # Apply the counting function to update the column
        self[npeaks_on] = self[ms_on].apply(count_peaks)







    
        
#################################################################
    def classfy_np(self, smiles_on='SMILES'):
        df = np_classfy_df(self, smiles_on=smiles_on)
        return df
    
    def drop_istd_peak(self, mz, rt, mz_window = 0.005, rt_window = 3):
        '''
        drop internal standard according to precursor mz and retentontime

        param:
            self, PeakFrame or pandas data frame object
            mz, precursor mz
            rt, retention time (min)
        return
            result data frame
        '''
        istd = self[(self['precursormz'] < mz + mz_window) & 
                    (self['precursormz'] > mz - mz_window) &
                    (self['retentiontime'] < rt + rt_window) &
                    (self['retentiontime'] > rt - rt_window)]
        return self.drop(istd.index)
    
    def drop_duplicated_ms(self, 
                           mz_on='precursormz',
                           MSMS_on='MSMS',
                           tol=(0.003, 0.005),
                           sim_thd={'bonanza':0.9, 'entropy':0.9, 'matched_ratio': 0.25},
                           keep_first_on = None,
                           ascending=False):
        '''
        drop duplicated MSMS

        param:
            mz_on, MSMS_on:  columns names of precursor mz and MSMS spectra
            tol:           tolerance for match
            smililarity:   similarity threshold for duplicates judgement. 
                            Based on the kernel density analysis of the metabolomics data from zebrafish,
                            the cutoff value for the similarity of MSMS matches should be set between 
                            0.85 and 0.92.
            keep_first_on: if None, retain the first one in the order of appearance.
                            if not None, sort by the column name specified in this parameter, 
                                and then keep the first one.
            device:   determine to use cpu or gpu edition for Match class

        return
            a data frame after deduplicates.
        '''
        scores_names = {'matched_count', 'matched_ratio', 'bonanza', 'simple_dot', 'modified_dot', 'entropy'}
        keys = set(sim_thd.keys())
        if not keys <= scores_names:
            raise ValueError(f'keys unacceptable: {keys-scores_names}.\n{scores_names}')
        
        if keep_first_on:
            df = self.sort_values(by=keep_first_on, ascending=ascending).copy().reset_index()
        else:
            df = self.copy()

        scores = df.match(mz_on=mz_on,
                          MSMS_on=MSMS_on,
                          tol=tol)
        
        condition = True
        for key in sim_thd:
            condition = condition & (scores[key] > sim_thd[key])

        idx_to_drop = scores.loc[condition, 'que_idx'].unique().tolist()
        return df.drop(index=idx_to_drop)

    def extract_ion_chrom(self, target_mz,
            intensity_on='intensity',
            precursormz_on='precursormz',
            ms1_error = 0.003,
            thd_intensity = 0.02):
        '''
        extract EIC of target mz
        param:
            thd_intensity, thd_intensity * intensity as the cut off for intensity'''
        cdt1 = self[precursormz_on] < (target_mz + ms1_error)
        cdt2 = self[precursormz_on] > (target_mz - ms1_error) 
        eic  = self.loc[cdt1 & cdt2]
        intensity_max = eic[intensity_on].max()
        return eic[eic[intensity_on] > thd_intensity * intensity_max]
    
    def enrich(self,
               que,                     # que, peak dataframe
               target_on,               # column name of enrich targets
               mz_on = 'precursormz',
               MSMS_on = 'MSMS',
               que_mz_on = None,
               que_msms_on = None,
               tol = (0.003, 0.005),
               sim_thd= {'bonanza':0.9, 'entropy':0.9},
               test_method = 'fisher',
               use_fdr = True): # return enrich data frame
        '''
        
             # (matched_count, matched_ratio, bonanza, simple_dot, modified_dot, entropy)
        '''      

        if not isinstance(que, self.__class__):
            raise TypeError(f'que is not {self.__class__.__name__} object!')
        
        scores_names = {'matched_count', 'matched_ratio', 'bonanza', 'simple_dot', 'modified_dot', 'entropy'}
        keys = set(sim_thd.keys())
        if not keys <= scores_names:
            raise ValueError(f'keys unacceptable: {keys-scores_names}.\n{scores_names}')
                
        # 按 'tcm_name' 列分组
        grouped = self.groupby(target_on)
        # 统计匹配数
        matches = []
        
        total = len(grouped)
        for i, (key, df) in enumerate(grouped, 1):
            # key 就是当前分组的 target_on 值
            print(f"\rProcessing {i}/{total} ({100.0 * i / total:.1f}%) - {target_on}={key}", end="", flush=True)
            scores = df.match(que,
                              mz_on=mz_on,
                              MSMS_on=MSMS_on,
                              que_mz_on=que_mz_on,
                              que_msms_on=que_msms_on,
                              tol=tol)

            condition = True
            for key in sim_thd:
                condition = condition & (scores[key] > sim_thd[key]) 
           
            n_match = scores.loc[condition, 'idx'].nunique()
            matches.append({target_on: key,
                            'n_match': n_match}) 

        matches = pd.DataFrame(matches)
        matches = matches[matches['n_match'] > 0].sort_values(by='n_match', ascending=False)
        matches.set_index(target_on, inplace=True)

        n_feature = self[target_on].value_counts()
        n_feature = n_feature.to_frame(name='n_feature')
        n_totlal_features = self.shape[0]
        counts = matches.join(n_feature, how='left')
        n_total_matches = int(counts['n_match'].sum())

        if counts.empty or counts.shape[0] == 0:
            warnings.warn("No matches. Returned an empty data frame.", UserWarning)
            return counts

        enr = enrich_df(counts[['n_feature', 'n_match']],
                        n_totlal_features,
                        n_total_matches,
                        use_fdr=use_fdr,
                        method=test_method)
        return enr  


    def ex_ms_array_3d(self, ms_on='MSMS', T=False, ms_len=0, npeaks_on='Num Peaks', dtype=np.float64):
        '''
        Convert MSMS column (containing multiple 2D arrays) to a 3D numpy array
        
        Parameters:
        -----------
        ms_on : str, default='MSMS'
            Column name containing MSMS data (each element should be a 2D numpy array)
        
        T : bool, default=False
            Indicates whether to transpose output or input format:
            - If input is not transposed (N, 2) and T=False: output (N_records, N, 2)
            - If input is transposed (2, N) and T=False: output (N_records, N, 2) - will untranspose
            - If T=True: output (N_records, 2, N)
            
        ms_len : int or None, default=None
            Specify fixed length for MSMS arrays (same logic as convert_ms_array):
            - None: keep original length (no padding)
            - 0: use maximum value in npeaks_on column
            - positive int: use this value as fixed length
                If smaller than max peaks count, raise warning and truncate data
            
            When ms_len is set, insufficient peaks are padded with [0, 0]
        
        npeaks_on : str, default='Num Peaks'
            Column name containing number of peaks, used when ms_len=0
        
        dtype : numpy dtype, default=np.float64
            Data type for the output array
        
        Returns:
        --------
        numpy.ndarray
            3D array with shape depending on T parameter:
            - If T=False: (N_records, N_peaks, 2) - standard format
            - If T=True: (N_records, 2, N_peaks) - transposed format
        
        Raises:
        -------
        ValueError
            - If ms_on column not found
            - If ms_len is negative or invalid
            - If npeaks_on column not found when ms_len=0
            - If MSMS arrays have inconsistent dimensions
        
        Warnings:
        ---------
        UserWarning
            - If ms_len is smaller than max peaks count (data will be truncated)
        
        Examples:
        ---------
        >>> pf = PeakFrame(...)
        
        # Variable length, standard format (N_peaks, 2)
        >>> arr_3d = pf.to_array_3d()
        >>> arr_3d.shape
        (5, 15, 2)
        
        # Variable length, transposed format (2, N_peaks)
        >>> arr_3d = pf.to_array_3d(T=True)
        >>> arr_3d.shape
        (5, 2, 15)
        
        # Fixed length to max peaks, standard format
        >>> arr_3d = pf.to_array_3d(ms_len=0)
        >>> arr_3d.shape
        (5, 61, 2)
        
        # Fixed length to max peaks, transposed format
        >>> arr_3d = pf.to_array_3d(ms_len=0, T=True)
        >>> arr_3d.shape
        (5, 2, 61)
        
        # Fixed length with specific value
        >>> arr_3d = pf.to_array_3d(ms_len=50, T=True)
        >>> arr_3d.shape
        (5, 2, 50)
        '''
        if ms_on not in self.columns:
            raise ValueError(f"Column '{ms_on}' not found in PeakFrame")
        if npeaks_on not in self.columns:
            raise ValueError(f"Column '{npeaks_on}' not found in PeakFrame")

        df = self.convert_ms_array(ms_on=ms_on,
                                   T=T,
                                   inplace=False,
                                   ms_len=ms_len,
                                   npeaks_on=npeaks_on)
        msms_array = df['MSMS'].values
        n_records = len(msms_array)
        arr_3d = np.stack([msms_array[i] for i in range(n_records)], axis=0)
        arr_3d = np.stack(msms_array, axis=0)
        
        return arr_3d
        
            
        
       



    
    # def find_precursor_type(self, target_mass, ionmode, mz_on):
    #     '''
    #     正负离子分开处理，找到的信号再合并
    #     该函数要适合self中同时具有有正负离子信息、只具有正离子或只具有负离子的情况
    #     find out precursor type according to the target compound mass (target_mass)
    #     param:
    #         target_mass, mass of target compound
    #         ionmode, postive or negative
    #     returns:
    #         mzfram containing matched results.
    #     '''
    #     from .precursortype import load_precursors
    #     from . import mz
    #     df = self.copy()
    #     df['Num Peaks'] = df['Num Peaks'].astype(int)
    #     df = df[df['Num Peaks'] > 0]
    #     df['precursortype'] = ''
    #     pcs = load_precursors(target_mass, ionmode)
    #     pcs = pcs[pcs['mz'] > 70]
    #     for idx in df.index:
    #         for j in pcs.index:
    #             if mz.match(df.loc[idx, mz_on], pcs.loc[j, 'mz']) == True:                    
    #                 df.loc[idx,'precursortype'] = pcs.loc[j, 'type']
    #                 break
    #     return df[df['precursortype'] != '']
    

    def flatten_msms_mz(self, intensity_tol=0, MSMS_on='MSMS', num_peaks_on = 'Num Peaks'):
        '''
        Obtain a flat array consisting of all mz values in ms2
        param:
            intensity_tol, intensity tolerance in MSMS
            MSMS_on, the column name of MSMS
            num_peaks_on, the name of column of "Num Peaks"
        '''
        MSMS = self.loc[self[num_peaks_on] > 0, MSMS_on]

        mz_values = [pair[0]            # 取第一列元素
                    for sub in MSMS       # sub = [] 或 [[x1,y1], [x2,y2], ...]
                    for pair in sub
                    if pair[1] > intensity_tol]   # pair = [x, y]

        return np.array(mz_values, dtype=float)

    
    def match(self,
              que=None,
              mz_on = 'precursormz',
              MSMS_on = 'MSMS',
              que_mz_on = None,
              que_msms_on = None,
              tol = (0.003, 0.005)):
        '''
        Calculate the MSMS similarity matrix between two peak frames (self and que).

        if que is None, match self

        return:
            similarity matrix or long table
        '''           

        if que_mz_on is None:
            que_mz_on = mz_on
        if que_msms_on is None:
            que_msms_on = MSMS_on

        if mz_on not in self.columns:
            raise ValueError(f'not found the column name {mz_on}')
        if MSMS_on not in self.columns:
            raise ValueError(f'not found the columns name {MSMS_on}')
        if que is not None:
            if not isinstance(que, self.__class__):
                raise TypeError(f'que is not {self.__class__.__name__} object')
            else:            
                if que_mz_on not in que.columns:
                    raise ValueError(f'not found the columns name {que_mz_on}')
                if que_msms_on not in que.columns:
                    raise ValueError(f'not found the columns name {que_msms_on}')
        
        self_msl = similarity.join_array(self[mz_on].values, self[MSMS_on].values)
        if que is not None:
            que_msl = similarity.join_array(que[que_mz_on].values, que[que_msms_on].values)
        else:
            que_msl = None
            
        scores= similarity.get_scores_batch(self_msl, que_msl, tol)
        # (matched_count, matched_ratio, bonanza, simple_dot, modified_dot, entropy) by batch
        score_names = ['matched_counts', 'matched_ratio', 'bonanza',
                       'simple_dot', 'modified_dot', 'entropy']
        score_dfs = [
            pd.DataFrame(scores[i]).stack()
            .rename_axis(['idx', 'que_idx']).reset_index(name=score_names[i])
            for i in range(6)
        ]
        #! 保持顺序一致，self.enrich函数的计算逻辑，依赖这个顺序

        # 合并所有数据框
        df = score_dfs[0]
        for s_df in score_dfs[1:]:
            df = df.merge(s_df, on=['idx', 'que_idx'])
        ## 位置索引转换为行索引
        df['idx'] = self.index[df['idx'].tolist()]
        if que is None:
            df['que_idx'] = self.index[df['que_idx'].tolist()]
            return df[df['idx'] < df['que_idx']]
        else:
            df['que_idx'] = que.index[df['que_idx'].tolist()]
            return df
        
    def match_by_chunk_to_csv(self,
                          save_to,
                          que=None,
                          mz_on = 'precursormz',
                          MSMS_on = 'MSMS',
                          que_mz_on = None,
                          que_msms_on = None,
                          tol = (0.003, 0.005),
                          similarity_cutoff_to_save=0.3,
                          chunk_size=10000):
        '''
        超大表格match运算
        '''
        warnings.warn("This function has not been tested and verified.", UserWarning)

        with open(save_to, 'w') as csv:
            csv.write('idx,que_idx,matched_counts,bonanza,cosine\n')
            csv.flush()

            if que is not None:
                for i in range(0, len(self), chunk_size):
                    chunk_self = self.iloc[i:i+chunk_size]
                    for j in range(0, len(que), chunk_size):
                        chunk_que = que.iloc[j:j+chunk_size]
                        scores_df = chunk_self.match(que=chunk_que,
                                                     mz_on = mz_on,
                                                     MSMS_on = MSMS_on,
                                                     que_mz_on = que_mz_on,
                                                     que_msms_on = que_msms_on,
                                                     tol =tol)
                        scores_df = scores_df[scores_df['bonanza'] > similarity_cutoff_to_save]
                        scores_df.to_csv(csv, mode='a+', index=False, header=False)
            else:
                chunks = []  
                for start in range(0, len(self), 3000):  
                    chunk = self.iloc[start:start+3000]  
                    chunks.append(chunk) 
   
                for i, chunk in enumerate(chunks[:-1]):
                    scores_df = chunk.match(mz_on = mz_on,
                                            MSMS_on = MSMS_on,
                                            tol =tol)
                    scores_df = scores_df[scores_df['bonanza'] > similarity_cutoff_to_save]
                    scores_df.to_csv(csv, mode='a+', index=False, header=False)

                    for next_chunk in chunks[i+1:]:
                        scores_df = chunk.match(que=next_chunk,
                                                mz_on = mz_on,
                                                MSMS_on = MSMS_on,
                                                tol =tol)
                        scores_df = scores_df[scores_df['bonanza'] > similarity_cutoff_to_save]
                        scores_df.to_csv(csv, mode='a+', index=False, header=False)

                scores_df = chunks[-1].match(mz_on = mz_on,
                                             MSMS_on = MSMS_on,
                                             tol =tol)
                scores_df = scores_df[scores_df['bonanza'] > similarity_cutoff_to_save]
                scores_df.to_csv(csv, mode='a+', index=False, header=False)

            csv.write('# finished.')      
     
    def match_counts(self,
                     que,                     # que, peak dataframe
                     target_on,               # 指定计算匹配数目的目标列
                     mz_on = 'precursormz',
                     MSMS_on = 'MSMS',
                     que_mz_on = None,
                     que_msms_on = None,
                     tol = (0.003, 0.005),
                     sim_thd= 0.9,
                     sim_type='bonanza'):       # similarity cut off
        '''

        return 
            df: 每个target的特征数和命中特征数
            int: 以及命中特征总数        
        '''   

        scores = self.match(que,
                            mz_on=mz_on,
                            MSMS_on=MSMS_on,
                            que_mz_on=que_mz_on,
                            que_msms_on=que_msms_on,
                            tol=tol)

        scores_matched = scores[scores[sim_type] > sim_thd]

        n_matched = self.loc[scores_matched['idx'].tolist(), target_on].value_counts()
        n_feature = self[target_on].value_counts()

        n_matched = n_matched.to_frame(name='n_matched')
        n_feature = n_feature.to_frame(name='n_feature')

        counts = n_feature.join(n_matched, how='left')
        counts['n_matched'] = counts['n_matched'].fillna(0).astype(int)

        counts = counts[counts['n_matched'] > 0]

        return counts.sort_values(by='n_matched', ascending=False)
    
    def norm_msms(self, MSMS_on='MSMS', inplace=False):
        '''
        Convert MSMS peak intensities to the percentage relative to
            the strongest peak (relative intensity)
        '''
        print('The code hasn’t been completed yet;')
        print('nothing is being executed.')
        pass  
   
    def round_msms(self, MSMS_on='MSMS', n = 5):
        '''
        Specify mz decimal places in MSMS
        n, 小数的保留位数
        '''
        self[MSMS_on] = self[MSMS_on].apply(lambda x:
                [[round(mz, n), i] for mz, i in x])        
    
    def standardize(self):
        '''
        Generate a standardized DataFrame:
            remove extra spaces( ), hyphens (-), or underscores (_) from column names;
            keep only standardized columns;
            refresh num peaks
        '''
        df = self.copy()
        if 'Num Peaks' in df.columns:
            df.drop(columns='Num Peaks', inplace=True)
        df.columns = [re.sub(r'[ _-]+', '', c).upper() for c in df.columns]
  
        df.update_num_peaks() # 会重新产生标准列名Num Peaks

        # 严格按 __FIELDS__ 顺序筛选已有列
        ordered_cols = [c for c in __FIELDS__ if c in df.columns] 
        ordered_cols.append('Num Peaks')     

        return df[ordered_cols]

    def to_msp_block(self,
                     MSMS_on: str = "MSMS",
                     npeaks_on: str='Num Peaks',
                     MSMS_sep: str = "\t") -> str:
        """
        将 DataFrame 转换为格式化文本，针对严格的 MSMS 数据结构（二维列表/数组，每行2列）。

        规则：
        - 非 MSMS 列：输出 "列名:值" 每列一行（NaN 输出为空值，即 '列名:'）。
        - MSMS 列：每个 [mz, intensity] 输出为 "mz{MSMS_sep}intensity " 并换行。
        - 记录之间不插入分隔线。

        参数：
        - df: pandas DataFrame
        - MSMS_on: MS/MS 列名（默认 "MSMS"）
        - MSMS_sep: MSMS 中 mz 与 intensity 的分隔符（默认 tab，即 "\\t"）

        返回：
        - str: 格式化后的文本
        """
        # 按当前 DataFrame 列顺序，排除 MSMS 列
        normal_cols = [c for c in self.columns if c not in (MSMS_on, npeaks_on) ]

        out_lines = []
        append = out_lines.append  # 局部绑定以加速循环

        for _, row in self.iterrows():
            # 1) 非 MSMS 列
            for c in normal_cols:
                append(f"{c}: { row[c]}")                   

            # 2) MSMS 列（严格二维列表，每行2列）
            if  MSMS_on in self.columns:
                msms_list = np.asarray(row[MSMS_on])  # [[mz, intensity], ...]
                if (msms_list is None) or (msms_list.size == 0):
                    append('Num Peaks: 0')
                else:
                    append(f'Num Peaks: {len(msms_list)}')
                    for mz, intensity in msms_list:
                        append(f"{mz}{MSMS_sep}{intensity}")
            append('') # 每条记录以空行分隔

        return "\n".join(out_lines)

    def to_msp(self,
               filename,
               standardized: bool = False,
               mode='w',
               MSMS_on: str = "MSMS",
               npeaks_on='Num Peaks',
               MSMS_sep: str = "\t",
               chunk_size: int = 5000,
               encoding='utf-8'):
        """
        将全量数据按 chunk_size 分块，逐块转换为 MSP 文本并写入文件。

        参数：
        - filename: 目标文件路径
        - mode: 文件写入模式（'w' 覆盖写入，'a' 追加写入等）
        - MSMS_on: MS/MS 列名（默认 "MSMS"）
        - MSMS_sep: MSMS 中 mz 与 intensity 的分隔符（默认 tab，即 "\\t"）
        - chunk_size: 每块行数（建议 100~5000 范围内，根据内存与速度平衡）
        - encoding: 文件编码（默认 'utf-8'）
        - standardized: 是否先对数据进行标准化（调用 self.standardize()）再导出

        依赖：
        - 需要类中提供 `to_msp_block` 方法，将指定行区间转换为 MSP 文本。
        """
        # 选择数据源：标准化或原始
        df_to_write = self.standardize() if standardized else self

        # 参数基本校验
        if chunk_size <= 0:
            raise ValueError("chunk_size must be a positive integer")
        
        if not isinstance(filename, str):
            raise ValueError(f'unknown file name: {filename}')

        # 计算总行数
        try:
            total = len(df_to_write.df)
        except AttributeError:
            total = len(df_to_write)

        if total == 0:
            # 空数据直接创建/清空文件后返回
            with open(filename, mode=mode, encoding=encoding) as f:
                pass
            return

        # 逐块写入
        with open(filename, mode=mode, encoding=encoding) as f:
            for start in range(0, total, chunk_size):
                end = min(start + chunk_size, total)
                # 注意：这里基于 df_to_write 的切片调用 to_msp_block
                block_text = df_to_write.iloc[start:end].to_msp_block(MSMS_on=MSMS_on,
                                                                      npeaks_on=npeaks_on,
                                                                      MSMS_sep=MSMS_sep)

                block_text = block_text.rstrip('\n') + '\n\n'
                f.write(block_text)

    def to_pickle(self, path, MSMS_on='MSMS', to_msms_str=False, *args, **kwargs):
        df = self.copy()
        if to_msms_str and (MSMS_on in df.columns):
            df[MSMS_on] = df[MSMS_on].apply(str)            
        return super().to_pickle(path, *args, **kwargs)
    
    def to_sqlite3(self, tbl_name, conn, if_exists='replace', index=False, MSMS_on='MSMS'):
        df = self.copy()
        if MSMS_on in df.columns:
            df[MSMS_on] = df[MSMS_on].apply(mz.to_str) # 仅当msms是np数组时才有效
        return df.to_sql(tbl_name, conn, if_exists=if_exists, index=index)
   
    def update_num_peaks(self):
        '''
        Update num peak field (number of fragments
        '''
        if 'MSMS' in self.columns:
            self['Num Peaks'] = self['MSMS'].apply(len)
        elif 'msms' in self.columns:
            self['Num Peaks'] = self['msms'].apply(len)
        else:
           self['Num Peaks'] = 0
           warnings.warn('Warning: no MSMS found!', UserWarning)

### readers
#-----------------------------------------------------------------------------------

def read_mgf(fpath,
             sep_ms2=' ', 
             ionmode:str = 'auto',
             encoding='utf-8', keep_raw_data=False):
    '''
    params:
        ionmode, can be 'pos', 'neg' or 'auto'. 
                    If auto, ionmode will be set automatically according to CHARG in records.
    '''
    data = []
    with open(fpath, encoding=encoding) as f:
        lines = f.readlines()
        for line in lines:
            if line.strip() == 'BEGIN IONS':
                item = {}
                item['MSMS'] = []
            elif line.strip() == 'END IONS':
                item['Num Peaks'] = len(item['MSMS'])
                data.append(item)
            elif '=' in line:
                keys = line.split('=')
                item[keys[0].strip()] = keys[1].strip()
            elif line[0].isdigit():
                mz, intensity = line.strip().split(sep_ms2)
                item['MSMS'].append([float(mz), float(intensity)])
    
    mgf = PeakFrame(data)
    
    if 'CHARGE' in mgf.columns:
        mgf['CHARGE'] = mgf['CHARGE'].fillna('')
    
    if 'RTINSECONDS' in mgf.columns:
        mgf['RTINSECONDS']   = mgf['RTINSECONDS'].astype(float)
        mgf['retentiontime'] = mgf['RTINSECONDS']/60
    
    # 处理 PEPMASS 列，支持只有 mz 或 mz+intensity 的情况
    if 'PEPMASS' in mgf.columns:
        # 分割 PEPMASS，expand=True 会返回所有分割的列
        pepmass_split = mgf['PEPMASS'].str.split(' ', expand=True, n=1)
        
        # 第一列总是存在的（mz值）
        mgf['precursormz'] = pd.to_numeric(pepmass_split[0], errors='coerce')
        
        # 第二列可能不存在（intensity值）
        if len(pepmass_split.columns) > 1:
            mgf['intensity'] = pd.to_numeric(pepmass_split[1], errors='coerce')
        else:
            # 如果没有intensity值，用NaN填充
            mgf['intensity'] = np.nan
    
    if ionmode == 'auto':
        if 'CHARGE' in mgf.columns:
            mgf['ionmode'] = mgf['CHARGE'].apply(lambda x: \
                'Positive' if x.endswith('+') \
                    else ('Negative' if x.endswith('-') else ''))
        else:
            mgf['ionmode'] = ''
    else:
        mgf['ionmode'] = ionmode
    
    mgf['comment'] = fpath

    return mgf.reset_index(drop=True)


def read_mona_msp(fpath,
                  extract_smiles = True,
                  sep_ms2=' '):
    '''
    read mona msp file, not suitable for other msp file.
    param:
        extract_smiles, extract smiles string from comment field.
    '''
    df = PeakFrame.read_msp(fpath, sep_ms2=sep_ms2)
    if extract_smiles:
        df['smiles'] = df['Comments'].str.extract('SMILES=(.*?)"')
    return df

def read_msd_msp(fname, use_relative_abundance=False, include_fname=True, **kwargs):
    '''
    read peak list (msp format) exported from MS-Dial version 5.2 or higher
    peak height and peak area are transferred into relative value to the base peak.
    use_relative_abundance: Whether to Use Relative Abundance
    return:
        a PeakFrame
    '''
    def _ex_intensity_(comment):
        '''
        从注释文本(msp的comment字段)中提取峰高和峰面积
        '''
        pkheight = re.search(r'PEAKHEIGHT=(\d+)', comment or '')
        pkarea   = re.search(r'PEAKAREA=(\d+)', comment or '')
        pkid     = re.search(r'PEAKID=(\d+)', comment or '')

        # 修改逻辑：确保即使匹配不到也不会返回 None
        pk_id          = int(pkid.group(1)) if pkid else -1
        pkheight_value = float(pkheight.group(1)) if pkheight else np.nan
        pkarea_value   = float(pkarea.group(1)) if pkarea else np.nan

        return pk_id, pkheight_value, pkarea_value

    df = read_msp(fname, **kwargs)

    # 解析 COMMENT 列
    df[['pkid', 'peak_height', 'peak_area']] = df['COMMENT'].apply(_ex_intensity_).apply(pd.Series)

    # 数据类型转换
    df['pkid'] = df['pkid'].astype(int)
    df[['peak_height', 'peak_area']] = df[['peak_height', 'peak_area']].astype(float)

    if use_relative_abundance:
        base_pk_height = df['peak_height'].max(skipna=True)
        base_pk_area  = df['peak_area'].max(skipna=True)
        if base_pk_height and not np.isnan(base_pk_height):
            df['peak_height'] = 100 * df['peak_height'] / base_pk_height
        if base_pk_area and not np.isnan(base_pk_area):
            df['peak_area'] = 100 * df['peak_area'] / base_pk_area

    if include_fname:
        df['file'] = str(fname)

    return df

def read_msd_msp_folder(folder, use_relative_abundance=False, **kwargs):
    '''
    Read all MSP files in the specified folder and merge them into a single peakFrame.
    '''
    from pathlib import Path
    folder_path = Path(folder)
    msp_files = [f for f in folder_path.iterdir() if f.is_file() and f.suffix.lower() == '.msp']
    
    if len(msp_files) == 0:
        raise ValueError(f'not msp file was found in {folder}')
    data = []
    for f in msp_files:
        print(f'reading:\t{f}', end='\r', flush=True)
        df = read_msd_msp(str(f), use_relative_abundance, include_fname=True,  **kwargs)
        if not df.empty:
            data.append(df)
    
    return pd.concat(data, ignore_index=True)


def read_msp(fpath,
             sep_ms2='\t',
             rename: dict=None,
             comment=None,
             to_float: set = {'PRECURSORMZ','RETENTIONTIME', 'INTENSITY', 'Num Peaks'},
             encoding='utf-8'):
    # 使用 pandas 读取文本文件
    msp = pd.read_table(fpath,
                        dtype=str,
                        sep='\r',
                        skip_blank_lines=False,
                        comment=comment,
                        header=None,
                        names=['txt'],
                        engine='c',
                        encoding=encoding)

    # 添加辅助列 'group_id' 标识文本块
    msp['group_id'] = (msp['txt'].isnull().cumsum())

    # 填充空值以便于后续处理
    msp.fillna('', inplace=True)

    # 合并同一文本块的行
    txt_blocks = msp.groupby('group_id')['txt'].apply('\n'.join).reset_index()

    # 过滤掉空的文本块
    txt_blocks = [block for block in txt_blocks['txt'] if block]

    # 如果需要重命名字段，替换文本块中的旧名称
    if rename is not None:
        for old_name, new_name in rename.items():
            txt_blocks = [block.replace(old_name, new_name) for block in txt_blocks]

    # 解析文本块为离子对象并创建 DataFrame
    ions = [PeakFrame._parse_msp_txt(block, sep_ms2=sep_ms2) for block in txt_blocks]
    df = PeakFrame(ions)

    # 将指定列转换为浮点数类型
    for col in to_float.intersection(df.columns):
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    if 'COMMENT' in df.columns and 'Comment' in df.columns:
        # 同时存在COMMENT和Coment两列的情况处理
        df['COMMENT'] = (
            df['COMMENT'].fillna('') +
            ((' | ' + df['Comment'].astype(str)).where(df['Comment'].notna(), ''))
        )
        df.drop(columns=['Comment'], inplace=True)

    return df


def read_pickle(fname, MSMS_on='MSMS', force_msms=False):
    '''
    read pickle file of PeakFrame
    param:
        fname, pickle file name
        MSMS_on, column name for MSMS
        force_msms, Whether to force the parsing conversion of MSMS, 
            pd.read_pickle seems to automatically convert list strings to arrays.

    由于模块名更改，在导入以前旧的模块名保存的pickle文件时，再度加载会找不到原模块名而报错。
    解决方案是加入模块别名：
    from mzpy import mzPandas as mpd
    sys.modules['mzpy.PeakFrame'] = mpd
    '''
    df = pd.read_pickle(fname)
    if force_msms and (MSMS_on in df.columns):
        df[MSMS_on] = df[MSMS_on].apply(ast.literal_eval)
    return PeakFrame(df)


def read_sql(sql, conn, MSMS_on='MSMS', parsing_msms=True):
    '''
    the method pd.read_sql does not automatically convert array strings to arrays by default, 
        so force_msms default is True (open).
    '''
    df = pd.read_sql(sql, conn)
    if parsing_msms and (MSMS_on in df.columns):
        df[MSMS_on] = df[MSMS_on].apply(ast.literal_eval)
    return PeakFrame(df) 
