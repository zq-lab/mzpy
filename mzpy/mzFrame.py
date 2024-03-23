#enconding: UTF-8
'''
msms database or data sheet processor

MSP:
    In every filed, it cannot contain newline character, such as "\r" or "\n".
        Or the msp text structure will be interrupted after being output.
        Therefore, in the mzFrame.to_msp method, "\r" and "\n" were checked and deleted firstly.
    FORMULA can be ''. But it can not be "nan" which can not be accpted by MS-Dial.
    MS-Dial does not accept single autom or ion, such as Na, N, S. 
        Thus items also be checked atom bumber before being exported in to_msp function.
'''

import ast
import numpy as np
import pandas as pd
import random
import re
import sqlite3
from tqdm import tqdm
from typing import List
from zipfile import ZipFile


def to_centroid(spectrum: np.ndarray,
                window_threshold_rate: float=0.33,
                mz_slice_width=0.1,
                n_peaks_threshold = 1) -> List[List[float]]:
    if len(spectrum) == 0:
        return []
    if not isinstance(spectrum, np.ndarray):
        spectrum = np.array(spectrum)
    
    uplift = spectrum[1:] > spectrum[:-1]
    if not uplift[:, 0].all():
        # 按mz大小排序
        spectrum = spectrum[np.argsort(spectrum[:, 0]), :]
    if len(spectrum) <= n_peaks_threshold:
        return spectrum
    
    # 峰检测的向量化操作
    uplift = uplift[:, 1]
    downlift = spectrum[1:, 1] < spectrum[:-1, 1]
    peaks_index: List[int] = np.where(uplift[:-1] & downlift[1:])[0] + 1    
    result: List[List[int]] = [None] * peaks_index.shape[0]
    
    for n, pidx in enumerate(peaks_index):
        # 从各峰中心开始，向两侧搜索数据点
        window_size: int = 1                                                        # 搜索的窗口大小
        center_mz, intensity_sum = spectrum[pidx]                                   # 该峰中心处的 mz, # 该峰中心处的 intensity (用于加权求 mz)
        weighted_mz: float = center_mz * intensity_sum                              # 用于加权求 mz 
        intensity_threshold: float = intensity_sum * window_threshold_rate          # intensity 阈值, 窗口搜索在窗口边界强度低于阈值时结束
        lp: np.ndarray = spectrum[pidx - 1]     # 窗口左边界的峰
        rp: np.ndarray = spectrum[pidx + 1]     # 窗口右边界的峰
        
        # 如果:
        # 窗口左边界的峰 intensiy 大于左边界左侧的峰 且
        # 窗口右边界的峰 intensiy 大于右边界右侧的峰 且
        # 窗口左边界与右边界的峰 intensity 均高于 intensity 阈值 且
        # 窗口左边界与右边界的峰 mz 与峰中心 mz 的偏差不超过 mz_slice_width
        # 则向左右扩展窗口        
        while pidx - window_size - 1 >= 0 and \
            pidx + window_size <= peaks_index.shape[0] - 2 and \
            uplift[pidx - window_size - 1] and downlift[pidx + window_size] and \
            (lp := spectrum[pidx - window_size - 1])[1] > intensity_threshold and \
            (rp := spectrum[pidx + window_size + 1])[1] > intensity_threshold and \
            abs(lp[0] - center_mz) < mz_slice_width and abs(rp[0] - center_mz) < mz_slice_width:           
            window_size += 1
            intensity_sum += lp[1] + rp[1]
            weighted_mz += lp[0] * lp[1] + rp[0] * rp[1]        
        # 计算加权 mz 后将该峰添加至结果中
        result[n] = [weighted_mz / intensity_sum, spectrum[pidx][1]]
    
    # 处理结果为空的特殊情况
    if not result:
        result: List[List[int]] = [spectrum[0], spectrum[-1]]        
    return result


      
class charge():
    pos = ('pos', 'positive', 'p', '+')
    neg = ('neg', 'negative', 'n', '-')

    @classmethod
    def standardize(cls, mode):
        if mode.lower() in cls.pos:
            return 'Positive'
        elif mode.lower() in cls.neg:
            return 'Negative'
        else:
            raise ValueError('not acceptable ionmode value')
        
    @classmethod
    def is_pos(cls, mode):
        return mode.lower() in cls.pos
    
    @classmethod
    def is_neg(cls, mode):
        return mode.lower() in cls.neg
     


class Precursor(): # 没大有存在的意义，应该并入mzFrame
    '''Standardize the keys of precursor ion to cooperate with MS-Dial
        - MS-Dial ignors case, such as ontology, when reading reference msp file.
        - MS-Dial ignors author and comment fields/
        - As a agreement, KEYS uses upper letters.
    '''

  
    __slots__ = ('name', 'precursormz', 'intensity', 'precursortype', 'formula',
                 'ontology', 'inchikey', 'smiles', 'retentiontime',
                 'ionmode', 'instrumenttype', 'instrument', 'collisionenergy',
                 'CCS', 'author', 'comment', 'msms')

    @classmethod
    def uniform_fields(cls):
        return cls.__slots__

    def __init__(self, data:dict = {}):
        for k in data:
            nk = re.sub(r'[^a-zA-Z]', '', k).lower()
            if nk in self.__slots__:
                self.__setattr__(nk, data[k])

    def __missing__(self, key):
        if isinstance(key, str):
            raise KeyError(key)
        return self[str(key)]

    def __contains__(self, key):
        return key in self.__slots__
    
    @classmethod
    def charge(cls, mode:str):
        if charge(mode) == 'Positive':
            return 'Positive'
        elif charge(mode) == 'Negative':
            return 'Negative'
    
    @classmethod
    def ms2_precision(cls, ms2:list):
        '''
        获取msms图谱的精度
        ms2:二维浮点数数组
        '''
        # 从第一列随机抽取三个元素，最多三个元素
        mz = random.sample([row[0] for row in ms2], min(len(ms2), 3))
        precision = []
        for data in mz:
            data_str = format(data)
            if '.' in data_str:
                precision.append(len(data_str)-data_str.index('.')-1)
            else:
                precision.append(0)
        return min(precision)

    def to_centroid(self,
                    window_threshold_rate: float=0.33,
                    mz_slice_width: float=0.1,
                    n_peaks_threshold:int = 1):
        self['msms'] = to_centroid(self['msms'],
                                   window_threshold_rate,
                                   mz_slice_width,
                                   n_peaks_threshold)

    def __str__(self, sep_ms2='\t'):
        txt = ''
        for field in self.__slots__:
            if field == 'msms':
                txt += f'Num Peaks: {len(self.msms)}\n'
                for data in self.msms:
                    txt += f'{data[0]}{sep_ms2}{data[1]}\n'
                txt += '\n'
            else:    
                txt += field.upper() + ': ' + str(getattr(self, field, '')) + '\n'
        return txt

    def to_dict(self):
        return {attr: getattr(self, attr, None) for attr in self.__slots__}

    @property
    def precision(self):
        if 'PRECURSORMZ' in self:
            return len(self['PRECURSORMZ'].split('.',1)[1])
        else:
            return -1
 
    def scale_msms(self, threshold = None):
        intensities = [i for mz, i in self.msms]
        base_intensity = max(intensities)
        self.msms = [[mz, round(100*i/base_intensity, 2)] for mz, i in self.msms]
        if threshold is not None:
            self.msms = [[mz, i] for mz, i in self.msms if i > threshold]



class mzFrame(pd.DataFrame):
    '''
    convention column name as Presursor.__slot__
    '''
    __count = 0

    @property
    def _constructor(self):
        return self.__class__
        # 使用pd.concat之后仍能保持子类类型
        # 使用pd.merge之后仍能保持子类类型
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        cols = {}
        for field in self.columns:
            k = re.sub(r'[^a-zA-Z]', '', field).lower()
            if k in Precursor.__slots__:
                cols[field] = k
        self.rename(columns=cols, inplace=True)
        # for col_name in (set(Precursor.__slots__) - set(self.columns)):
        #     self[col_name] = None

    ### internal parse    
    @classmethod
    def _parse_msp_txt(cls, txt, sep_ms2= '\t') -> dict:
        data = {}
        data['msms'] = []
        lines = txt.strip().split('\n')
        for line in lines:
            if ':' in line:
                d = line.split(':', )
                data[d[0]] = d[1].strip() 
            else:
                d = line.strip().split(sep_ms2, 1)
                data['msms'].append([float(d[0]), float(d[1])])
        return data
    
    @classmethod
    def _parse_cfmid_txt(cls, txt):
        ion = {}
        msms = []
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
            elif line.strip() == '': # 放在前面可以防止下一个判断的游标溢出
                ion['collisionenergy'] = '10, 20 40 V'                
                break
            elif line[0].isdigit():
                mz, intensity = line.split(' ')[0:2]
                msms.append([float(mz), float(intensity)])
        ion['msms'] = msms
        return ion
    
    ### reader method
    @classmethod
    def read_msp(cls,fpath, sep_ms2='\t',
                 rename:dict = None, 
                 comment=None,
                 to_float:list = {'PRECURSORMZ','RETENTIONTIME', 'INTENSITY'},
                 encoding='utf-8'):
        msp = pd.read_table(fpath,
                            dtype = str,
                            sep = '\r', # '\n'会报错，改用'/n'即可
                            skip_blank_lines = False, #! 避免跳过空行
                            comment = comment,
                            header = None,
                            names = ['txt'],
                            engine='c',
                            encoding = encoding)
        '''
        raname fields in msp txt, key as old name, value as new name
        '''
        cls.__count = 0
        def count_(x):
            if pd.isnull(x):
                cls.__count += 1
                return cls.__count-1
            else:
                return cls.__count
        msp['id'] = msp['txt'].apply(count_)
        cls.__count = 0
        msp.fillna('',inplace=True) #
        txt = msp.groupby('id')['txt'].apply('\n'.join).reset_index()
        txt = [x for x in txt['txt'] if x != '']
        if rename is not None:
            for k in rename:
                txt = [x.replace(k, rename[k]) for x in txt]
        ions = [cls._parse_msp_txt(x, sep_ms2=sep_ms2) for x in txt]
        df =  cls(ions)
        if 'precursormz' in df.columns:
            df['precursormz'] = df['precursormz'].astype(float)
        if 'retentiontime' in df.columns:
            df['retentiontime'] = df['retentiontime'].astype(float)
        if 'intensity' in df.columns:
            df['intensity'] = df['intensity'].astype(float)
        df['comment'] = fpath
        return df
    
    @classmethod
    def read_mgf(cls, fpath, sep_ms2=' ', 
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
                    item['msms'] = []
                elif line.strip() == 'END IONS':
                    item['Num Peaks'] = len(item['msms'])
                    data.append(item)
                elif '=' in line:
                    keys = line.split('=')
                    item[keys[0].strip()] = keys[1].strip()
                elif line[0].isdigit():
                    mz, intensity = line.strip().split(sep_ms2)
                    item['msms'].append([float(mz), float(intensity)])
        mgf = cls(data)
        if 'CHARGE' in mgf.columns:
            mgf['CHARGE'] = mgf['CHARGE'].fillna('')
        if 'RTINSECONDS' in mgf.columns:
            mgf['RTINSECONDS']   = mgf['RTINSECONDS'].astype(float)
            mgf['retentiontime'] = mgf['RTINSECONDS']/60
        mgf[['precursormz', 'intensity']] = mgf['PEPMASS'].str.split(' ', expand=True)
        mgf['precursormz'] = mgf['precursormz'].astype(float)
        mgf['intensity']   = mgf['intensity'].astype(float)
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

    @classmethod
    def read_mona_msp(cls, fpath,
                        extract_smiles = True,
                        sep_ms2=' ',
                        rename_columns = True):
        '''read mona msp file, not suitable for other msp file.'''
        df = cls. read_msp(fpath, sep_ms2=sep_ms2)
        if extract_smiles:
            df['smiles'] = df['Comments'].str.extract('SMILES=(.*?)"')
        return df
    
    @classmethod
    def read_cfmid_zip(cls, filename, unify_columns:bool = False, test = False, test_n = 5):
        data = []
        with ZipFile(filename, 'r') as zf:
            if test:
                namelist = random.sample(zf.namelist(), test_n)
            else:
                namelist = zf.namelist()
            for f in tqdm(namelist, desc='loading...', ncols=100):
                with zf.open(f, 'r') as file:
                    txt = file.read().decode()
                    ion = cls._parse_cfmid_txt(txt)
                    data.append(ion)
        print(f'Finished loading {filename}.')
        df = cls(data)
        return df

    @classmethod
    def read_sqlite3(cls, tbl_name, db_file, *args, **kwargs):
        # df = pd.read_sql(f'SELECT * FROM {tbl_name}', con, *args, **kwargs)
        # if 'msms' in df:
        #     df['msms'] = df['msms'].apply(ast.literal_eval)
        # return cls(df)
        with sqlite3.connect(db_file) as conn:
            df = pd.read_sql(f'SELECT * FROM {tbl_name}', conn, *args, **kwargs)
            if 'msms' in df:
                df['msms'] = df['msms'].apply(ast.literal_eval)
        return cls(df)

    
    @classmethod
    def read_pickle(cls, f, **kwargs):
        df = pd.read_pickle(f, **kwargs)
        return cls(df)
    
    def update_from_dict(self, data:dict):
        for key in data:
            self[key] = data[key] 
    
    ### gets and finds    
    def find_cols(self, name:str, return_name:bool = True):
        '''Find columns that match the name, 
            ignoring case and non-alphabetic characters.
        param:
            return_name: if true, return the matched column name(s)
                if false, returen the column objedct(s)
        '''
        name = ''.join(filter(str.isalpha, name.lower()))
        result = []
        for it in self.columns:               
            if name == ''.join(filter(str.isalpha, it.lower())):
                result.append(it)
        if return_name:
            return result
        else:
            return self[result]

    def get_top_peaks(self, top_n = 10):
        '''
        pick out peak (top) mz value 
        top_n, n top high intensity
        '''
        if ('precursormz' in self.columns) and ('intensity' in self.columns):
            peaks = to_centroid(self[['precursormz', 'intensity']].values)
            peaks = sorted(peaks, key=lambda x: x[1], reverse=True)
            seletect_mz = [it[0] for it in peaks[:top_n]]
            return self[self['precursormz'].isin(seletect_mz)]
        else:
            raise ValueError('The mz frame must contains "precursormz" and "intensity" columns')
     
    def select_ion(self, mz, on='precursormz'):
        error = (self[on]-mz)/mz
        return self[error < 5E-6]

    def select_pos(self, ionmode_on = 'ionmode'):
        selected = self[ionmode_on].apply(charge.is_pos)
        return self.loc[selected,]
    
    def select_neg(self, ionmode_on = 'ionmode'):
        selected = self[ionmode_on].apply(charge.is_neg)
        return self.loc[selected,]
    
    def filter_msms(self, thd_int:float = 0):
        self['msms'] = self['msms'].apply(lambda x:
                [[mz, intensity] for mz, intensity in x if intensity > thd_int])

    def centroid_msms(self):
        self['msms'] = self['msms'].apply(lambda x: to_centroid(x))
    
    def scale_msms(self, threshold = 0):
        '''
        Convert msms peak intensities to the percentage relative to
            the strongest peak (relative intensity)
        '''
        def percentage(x_y_data, threshold = 0):
            y_max = max([y for _, y in x_y_data])
            xy = [[x, round(100*y/y_max, 2)] for x, y in x_y_data]
            return [[x, y] for x, y in xy if y > threshold]                        
        self['msms'] = self['msms'].apply(lambda xy: percentage(xy, threshold))

    def unify(self):
        cols = {}
        for field in self.columns:
            k = re.sub(r'[^a-zA-Z]', '', field).lower()
            if k in Precursor.__slots__:
                cols[field] = k
        if 'Num Peaks' in self.columns:
            return self.rename(columns=cols)[list(cols.values()) + ['Num Peaks']]
        else:
            return self.rename(columns=cols)[list(cols.values())]

    def round_msms(self, n = 5):
        '''
        Specify mz decimal places
        n, 小数的保留位数
        '''
        self['msms'] = self['msms'].apply(lambda x:
                [[round(mz, n), i] for mz, i in x])


    ### writer and exporting
    def to_msp(self, filename, mode='w',
               check_smiles = False,
                chunk_size:int = 100,
                encoding='utf-8'):
        print('check data ...rows with null smiles or formula are dropped.')
        df = self.replace('\r', '', regex=True).replace('\n', '', regex=True)
        df = df.dropna(subset=['smiles', 'formula'])

        i = 0
        f = open(filename, mode=mode, encoding=encoding)
        for idx in tqdm(df.index,desc='writing...', ncols=100):
            i += 1
            ion = Precursor(df.loc[idx].to_dict())
            f.writelines(str(ion))
            if i%chunk_size == 0:
                f.close()
                f = open(filename, mode='a', encoding=encoding)
                i = 0
        f.close()

    def to_sqlite3(self, *args, **kwargs):
        df = self.unify()
        df['msms'] = df['msms'].apply(str)
        return df.to_sql(*args, **kwargs)
    

    ### plot chromatography
    def plot(self, x = 'retentiontime', y = 'intensity',
         legend= False, linewidth = 0.5,
         *args, **kwargs):
        return super().plot(x = x,
                            y = y,
                            legend = legend,
                            linewidth = linewidth,
                            *args,
                            **kwargs)


    def eic(self, target_mz, ms1_error = 0.003, thd_intensity = 0.02):
        '''
        extract EIC of target mz
        param:
            thd_intensity, thd_intensity * intensity as the cut off for intensity'''
        cdt1 = self['precursormz'] < (target_mz + ms1_error)
        cdt2 = self['precursormz'] > (target_mz - ms1_error) 
        eic  = self.loc[cdt1 & cdt2]
        intensity_max = eic['intensity'].max()
        return eic[eic['intensity'] > thd_intensity * intensity_max]

