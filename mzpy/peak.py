#enconding: UTF-8
'''
msms database or data sheet processor

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
from tqdm import tqdm

from . import mz
from .msl import MSList as MSL
from .stat import enrich_df


### basic function for PeakFrame
def concat(mzframe_list, msms_on='msms', ignore_index=False):
    for df in mzframe_list:
        df[msms_on] = df[msms_on].apply(mz.to_str)
    rst = pd.concat(mzframe_list, ignore_index=ignore_index)
    rst[msms_on] = rst[msms_on].apply(ast.literal_eval)   
    
    return rst

### operate functions for msms
def _ex_intensity_(comment):
    '''
    从注释文本(msp的comment字段)中提取峰高和峰面积
    '''
    pkheight = re.search(r'PEAKHEIGHT=(\d+)', comment)  
    pkarea = re.search(r'PEAKAREA=(\d+)', comment) 
    pkid =  re.search(r'PEAKID=(\d+)', comment) 

    # 获取提取的数值  
    pkheight_value = pkheight.group(1) if pkheight else None  
    pkarea_value   = pkarea.group(1) if pkarea else None
    pk_id          = pkid.group(1) if pkarea else -1
    return pk_id, float(pkheight_value), float(pkarea_value)


class Ion():
    '''Standardize the keys of precursor ion to cooperate with MS-Dial
        - MS-Dial ignors case, such as ontology, when reading reference msp file.
        - MS-Dial ignors author and comment fields/
        - As a agreement, KEYS uses upper letters.
    '''  
    __slots__ = ('name',       'precursormz',     'precursortype', 'ionmode',  'retentiontime',
                 'formula',    'ontology',        'smiles',        'inchikey', 'instrumenttype',
                 'instrument', 'collisionenergy', 'CCS',           'comment',  'msms')
    # 15 parameters in total
    # num peaks: is a derived value and does not appear a attributes.
    
    def __init__(self,
                 name:str = None,   precursormz:float=None,   precursortype:str=None, ionmode:str=None,
                 retentiontime:float=None, formula:str=None,  ontology:str=None,
                 smiles:str=None,   inchikey:str=None,        instrumenttype:str=None,
                 instrument:str=None,collisionenergy:str=None, CCS:float=None,
                 comment:str=None,   msms:list=None, *args, **kwargs):
        '''
        omit other parameters.
        '''
        self.name            = name
        self.precursormz     = precursormz
        self.precursortype   = precursortype
        self.retentiontime   = retentiontime
        self.formula         = formula
        self.ontology        = ontology
        self.smiles          = smiles
        self.inchikey        = inchikey
        self.instrumenttype  = instrumenttype
        self.instrument      = instrument
        self.collisionenergy = collisionenergy
        self.CCS             = CCS
        self.comment         = comment        
        if str(ionmode).lower() in ('+', 'p', 'pos', 'positive'):
            self.ionmode = 'Positive'
        elif str(ionmode).lower() in ('-', 'n', 'neg', 'negative'):
            self.ionmode = 'Negative'

        self.msms = mz.normalize(msms)   

    def __contains__(self, key):
        return key in self.__slots__
    
    def __missing__(self, key):
        if isinstance(key, str):
            raise KeyError(key)
        return self[str(key)]
    
    def __repr__(self) -> str:
        return 'Ion object \n---------------\n' + self.__str__()
    
    def __str__(self, sep_ms2='\t'):
        '''
        generate msp text block
        '''
        txt = ''
        for field in self.__slots__:
            if field == 'msms':
                if self.msms is None:
                    txt += 'Num Peaks: 0'
                else:
                    txt += f'Num Peaks: {len(self.msms)}\n'
                    for data in self.msms:
                        txt += f'{data[0]}{sep_ms2}{data[1]}\n'
                    txt += '\n'
            else:    
                txt += field.upper() + ': ' + str(getattr(self, field, '')) + '\n'
        return txt

    @classmethod
    def create_by_dict(cls, data:dict = {}):
        ion = cls()
        for k in data:
            nk = re.sub(r'[^a-zA-Z]', '', k).lower()
            if nk in ion.__slots__:
                ion.__setattr__(nk, data[k])  
    def centroid(self,
                    window_threshold_rate: float=0.33,
                    mz_slice_width: float=0.1,
                    n_peaks_threshold:int = 1):
        return mz.centroid(self.msms,
                            window_threshold_rate,
                            mz_slice_width,
                            n_peaks_threshold)  
    
    def drop_ms2_below(self, drop_min=0):
        self.msms = self.msms[self.msms[:, 1] >= drop_min]

    def match_mz(self, mz, tol=0.003, tol_rel=5, mode='abs'):
        '''
        Determine if two mz values (mz1 and mz2) match.
        param:
            mz1, mz2, mz values to compare.
            tol, tolerance. If mode is set as 'both', both tol and tol_rel are required.
        return:
            True or False
        '''
        return mz.match_mz(self.precursormz, mz, tol=tol, tol_rel=tol_rel, mode=mode)

    @property
    def npk(self):
        '''
        Num Peaks
        '''
        return len(self.ms2)    

    @property
    def precision(self):
        if 'PRECURSORMZ' in self:
            return len(str(self['PRECURSORMZ']).split('.',1)[1])
        else:
            return -1 
     
    def to_dict(self):
        return {slot: getattr(self, slot) for slot in self.__slots__}    

    def to_msp(self, sep_ms2='\t'):
        return self.__str__(sep_ms2=sep_ms2)
            
 

class PeakFrame(pd.DataFrame):
    '''
    convention column name as Presursor.__slot__
    '''

    @property
    def _constructor(self):
        return self.__class__
        # 使用pd.concat之后仍能保持子类类型
        # 使用pd.merge之后仍能保持子类类型
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # if transfer_col_names: # 试图添加这个参数，但是报错：__init__(self, *args, **kwargs)
        #     '''trnansfer colnames automanually'''
        #     cols = {}
        #     for field in self.columns:
        #         k = re.sub(r'[^a-zA-Z]', '', field).lower()
        #         if k in Ion.__slots__:
        #             cols[field] = k
        #     self.rename(columns=cols, inplace=True)
 
    @classmethod
    def _parse_msp_txt(cls, txt, sep_ms2= '\t') -> dict:
        data = {}
        data['msms'] = []
        lines = txt.strip().split('\n')
        for line in lines:
            if ':' in line:
                d = line.split(':', 1)
                data[d[0]] = d[1].strip() 
            elif line.strip(): # 只处理非空行非空行，但不含:
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
    
    def centroid_msms(self):
        self['msms'] = self['msms'].apply(lambda x: mz.centroid(x))


    ### plot chromatography
    def chrom(self, x = 'retentiontime', y = 'intensity',
              legend= False, linewidth = 0.5,
              *args, **kwargs):
        return super().plot(x = x,
                            y = y,
                            legend = legend,
                            linewidth = linewidth,
                            *args,
                            **kwargs)  
  
    def drop_istd(self, mz, rt, mz_window = 0.005, rt_window = 3):
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
                    (self['retentiontime'] > rt + rt_window)]
        return self.drop(istd.index)
    
    def drop_duplicated_ms(self, 
                           mz_on='precursormz',
                           msms_on='msms',
                           tol=(0.003, 0.005),
                           precursormz_compared=True,
                           similarity=0.99,
                           keep_first_on = None,
                           n_thread=1):
        '''
        drop duplicated msms

        param:
            mz_on, msms_on:  columns names of precursor mz and MSMS spectra
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
        if keep_first_on:
            df = self.sort_values(by=keep_first_on).copy().reset_index()
        else:
            df = self.copy().reset_index()

        msl = MSL(df[mz_on], df[msms_on])
        
        score = msl.compute_similarity_self(tol, precursormz_compared, n_thread=n_thread)
        
        upper_triangle = np.triu(score, k=1)
        _, cols = np.where(upper_triangle >= similarity)

        return df.drop(index=df.index[cols])

    def eic(self, target_mz,
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
               msms_on = 'msms',
               que_mz_on = None,
               que_msms_on = None,
               tol = (0.003, 0.005),
               device = 'cpu',
               similarity = 0.85,       # similarity cut off
               test_method = 'fisher',
               fdr = True

        ): # return enrich data frame
        '''
        Match based on the self ions and que ions,
            and enrich by rows according to the matches.
        '''      
   
        enr, n_total_hit = self.match_counts(que=que,
                                             target_on=target_on,
                                             mz_on=mz_on,
                                             msms_on=msms_on,
                                             que_mz_on=que_mz_on,
                                             que_msms_on=que_msms_on,
                                             tol=tol,
                                             device=device,
                                             similarity=similarity)
        
        enr = enrich_df(enr,
                        self.shape[0],
                        n_total_hit,
                        fdr=fdr,
                        method=test_method)
        
        return enr
    
    def find_precursor_type(self, target_mass, ionmode):
        '''
        find out precursor type according to the target compound mass (target_mass)
        param:
            target_mass, mass of target compound
            ionmode, postive or negative
        returns:
            mzfram containing matched results.
        '''
        from .precursorType import load_precursors
        df = self.copy()
        df['Num Peaks'] = df['Num Peaks'].astype(int)
        df = df[df['Num Peaks'] > 0]
        df['precursortype'] = ''
        pcs = load_precursors(target_mass, ionmode)
        pcs = pcs[pcs['mz'] > 70]
        for idx in df.index:
            mz = df.loc[idx, 'precursormz']
            for j in pcs.index:
                if mz.match_mz(mz, pcs.loc[j, 'mz']) == True:                    
                    df.loc[idx,'precursortype'] = pcs.loc[j, 'type']
                    break
        return df[df['precursortype'] != '']
    
    def match(self,
              que,
              mz_on = 'precursormz',
              msms_on = 'msms',
              que_mz_on = None,
              que_msms_on = None,
              tol = (0.003, 0.005),
              n_thread=1,
              return_matrix = False):
        '''
        Calculate the MSMS similarity matrix between two peak frames (self and que).

        return:
            similarity matrix or long table
        '''           

        if que_mz_on is None:
            que_mz_on = mz_on
        if que_msms_on is None:
            que_msms_on = msms_on

        if mz_on not in self.columns:
            raise ValueError(f'not found the column name {mz_on}')
        if msms_on not in self.columns:
            raise ValueError(f'not found the columns name {msms_on}')
        if que_mz_on not in que.columns:
            raise ValueError(f'not found the columns name {que_mz_on}')
        if que_msms_on not in que.columns:
            raise ValueError(f'not found the columns name {que_msms_on}')

        self_msl = MSL(self[mz_on], self[msms_on])
        que_msl  = MSL(que[que_mz_on], que[que_msms_on])
        
        scores = self_msl.compute_similarity(que_msl, tol=tol, n_thread=n_thread)
        if return_matrix:
            return scores # 此返回值中，不含有self和que的行索引
        else:
            # 获取行索引和列索引  
            row_indices, col_indices = np.indices(scores.shape)  

            # 创建一个DataFrame，使用 df1 和 df2 的索引  
            long_table = pd.DataFrame({  
                'index'     : self.index[row_indices.ravel()],   # df1 的行索引值(不是位置索引) 
                'que_index' :  que.index[col_indices.ravel()],   # df2 的列索引值  
                'similarity': scores.ravel()                             # 对应的值  
            })         
            return long_table
        
    def match_counts(self,
                     que,                     # que, peak dataframe
                     target_on,               # column name of enrich targets
                     mz_on = 'precursormz',
                     msms_on = 'msms',
                     que_mz_on = None,
                     que_msms_on = None,
                     tol = (0.003, 0.005),
                     device = 'cpu',
                     similarity = 0.85):       # similarity cut off
        '''
        计算self和que的中每对离子的相似度

        return 
            df: 每个target的特征数和命中特征数
            int: 以及命中特征总数        
        '''        
        base = self.set_index(target_on, drop=True)
        scores = base.match(que=que,
                            mz_on=mz_on,
                            msms_on=msms_on,
                            que_mz_on=que_mz_on,
                            que_msms_on=que_msms_on,
                            tol=tol,
                            device=device)
        hit_df = scores[scores['similarity'] > similarity]
        # 每个target的命中数统计
        hit = hit_df['index'].value_counts().to_frame(name='num_hit_features')
        # 每个target的features数目
        features = base.index.value_counts().to_frame(name='num_features')
        # counts汇总表
        counts = features.join(hit, how='left')
        counts['num_hit_features'] = counts['num_hit_features'].fillna(0).astype(int)
        
        # return counts[counts['num_hit_features'] > 0], hit_df['que_index'].nunique() 
        return counts, hit_df['que_index'].nunique()
        # 求算匹配离子总数时，que的一个离子有可能匹配到self的两个离子上
        # 因此不使用 counts['num_hit_features']计算

    def match_precursor_mz(self, mz, mz_on='precursormz', tol=0.003, tol_rel=5, mode='abs'):
        '''
        return precursor mz matched result
        '''
        condition = self[mz_on].apply(lambda x: mz.match_mz(x, mz, tol=tol, tol_rel=tol_rel, mode=mode))
        return self[condition]
    
    def match_self(self,
                   mz_on,
                   msms_on,
                   tol=(0.003, 0.005),
                   device = 'cpu',
                   return_matrix = False):
        if device == 'gpu':
            from .mslcp import MSList_cp as MSL
        else:
            from .msl import MSList as MSL

        if mz_on not in self.columns:
            raise ValueError(f'not found the column name {mz_on}')
        if msms_on not in self.columns:
            raise ValueError(f'not found the columns name {msms_on}')

        self_msl = MSL(self[mz_on], self[msms_on])
       
        scores = self_msl.compute_similarity_self(tol=tol)
        if return_matrix:
            return scores # 此返回值中，不含有self和que的行索引
        else:
            # 获取行索引和列索引  
            row_indices, col_indices = np.indices(scores.shape)  

            # 创建一个DataFrame，使用 df1 和 df2 的索引  
            long_table = pd.DataFrame({  
                'index'     : self.index[row_indices.ravel()],   # df1 的行索引值(不是位置索引) 
                'que_index' : self.index[col_indices.ravel()],   # df2 的列索引值  
                'similarity': scores.ravel()                             # 对应的值  
            })         
            return long_table
    
    def norm_msms(self, precursormz_on='precursormz', msms_on='msms',
                  Drop_unbroken_precursor=False,
                  tol=0.003, tol_rel=5, mode='abs'):
        '''
        Convert msms peak intensities to the percentage relative to
            the strongest peak (relative intensity)
        '''
        df = self.copy()
        if 'Num Peaks' not in df.columns:
            df['Num Peaks'] = df[msms_on].apply(len)
        df = df[df['Num Peaks'] > 0]
        
        idx_to_drop = []

        if Drop_unbroken_precursor:
            # 舍弃没有裂解的母离子
            for idx in df.index:
                if df.loc[idx, 'Num Peaks'] == 1:
                    if mz.match_mz(df.loc[idx, precursormz_on], df.loc[idx, msms_on][0][0],
                               tol=tol,
                               tol_rel=tol_rel,
                               mode=mode):
                        idx_to_drop.append(idx)                
            df = df.drop(idx_to_drop)
        
        df[msms_on] = df[msms_on].apply(lambda msms: mz.normalize(msms))
        return df    
   
    def round_msms(self, n = 5):
        '''
        Specify mz decimal places
        n, 小数的保留位数
        '''
        self['msms'] = self['msms'].apply(lambda x:
                [[round(mz, n), i] for mz, i in x])

    @property
    def top_peaks(self, top_n = 100, mz_on='precursormz', intensity_on='intensity'):
        '''
        pick out peak (top) mz value 
        top_n, n top high intensity
        '''
        top_n = min(top_n, self.shape[0])
        peaks = mz.centroid(self[[mz_on, intensity_on]].values)
        peaks = sorted(peaks, key=lambda x: x[1], reverse=True)
        seletect_mz = [it[0] for it in peaks[:top_n]]
        return self[self[mz_on].isin(seletect_mz)] 

    def to_msp_block(self):
        txt = ''
        for idx in self.index:
            ion = Ion(**self.loc[idx].to_dict())
            txt += str(ion)
        return(txt.strip()+'\n\n') 

    ### writer and exporting
    def to_msp(self, filename, mode='w',
                chunk_size:int = 500,
                encoding='utf-8'):
        print('check data ...rows with null smiles or formula are dropped.')
        df = self.replace('\r', '', regex=True).replace('\n', '', regex=True)
        # df = df.dropna(subset=['smiles', 'formula'])
            ## 会报错：expected string or bytes-like object, got 'int'

        df = df.reset_index()
        f = open(filename, mode=mode, encoding=encoding)
        for idx in tqdm(df.index, desc='writing...', ncols=100):
            ion = Ion(**df.loc[idx].to_dict())
            f.writelines(str(ion))
            if idx%chunk_size == 0:
                f.flush()
        f.close()

    def to_pickle(self, path, msms_on='msms', to_msms_str=False, *args, **kwargs):
        df = self.copy()
        if to_msms_str and (msms_on in df.columns):
            df[msms_on] = df[msms_on].apply(str)            
        return super().to_pickle(path, *args, **kwargs)
    
    def to_sqlite3(self, tbl_name, conn, if_exists='replace', index=False, msms_on='msms'):
        df = self.copy()
        if msms_on in df.columns:
            df[msms_on] = df[msms_on].apply(mz.to_str) # 仅当msms是np数组时才有效
        return df.to_sql(tbl_name, conn, if_exists=if_exists, index=index)


    def uniform(self):
        cols = {}
        for field in self.columns:
            k = re.sub(r'[^a-zA-Z]', '', field).lower()
            if k in Ion.__slots__:
                cols[field] = k
        if 'Num Peaks' in self.columns:
            return self.rename(columns=cols)[list(cols.values()) + ['Num Peaks']]
        else:
            return self.rename(columns=cols)[list(cols.values())]



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
    mgf = PeakFrame(data)
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


def read_mona_msp(fpath,
                  extract_smiles = True,
                  sep_ms2=' '):
    '''
    read mona msp file, not suitable for other msp file.
    param:
        extract_smiles, extract smiles string from comment field.
    '''
    df = PeakFrame. read_msp(fpath, sep_ms2=sep_ms2)
    if extract_smiles:
        df['smiles'] = df['Comments'].str.extract('SMILES=(.*?)"')
    return df


def read_msd_ali(fname, washed=False):
    from .metab import Metab
    return Metab.read_msd_alignment(fname, washed=washed)


def read_msd_msp(fname, rel_abd=False, **kwargs):
    '''
    read peak list (msp format) exported from MS-Dial version 5.2 or higher
    peak height and peak area are transfered into relative value to the base peak.
    rel_abd: Whether to Use Relative Abundance
    return:
        a PeakFrame
    '''
    df = read_msp(fname, ** kwargs)
    df[['pkid', 'peak_height', 'peak_area']] = df['COMMENT'].apply(_ex_intensity_).apply(pd.Series)

    if rel_abd:
        base_pk_heght = df['peak_height'].max()
        base_pk_area  = df['peak_area'].max()
        df['peak_height'] = 100*df['peak_height']/base_pk_heght
        df['peak_area']   = 100*df['peak_area']/base_pk_area

    return df


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

    return df


def read_pickle(fname, msms_on='msms', force_msms=False):
    '''
    read pickle file of PeakFrame
    param:
        fname, pickle file name
        msms_on, column name for MSMS
        force_msms, Whether to force the parsing conversion of MSMS, 
            pd.read_pickle seems to automatically convert list strings to arrays.

    由于模块名更改，在导入以前旧的模块名保存的pickle文件时，再度加载会找不到原模块名而报错。
    解决方案是加入模块别名：
    from mzpy import mzPandas as mpd
    sys.modules['mzpy.PeakFrame'] = mpd
    '''
    df = pd.read_pickle(fname)
    if force_msms and (msms_on in df.columns):
        df[msms_on] = df[msms_on].apply(ast.literal_eval)
    return PeakFrame(df)


def read_sql(sql, conn, msms_on='msms', force_msms=True):
    '''
    the method pd.read_sql does not automatically convert array strings to arrays by default, 
        so force_msms default is True (open).
    '''
    df = pd.read_sql(sql, conn)
    if force_msms and (msms_on in df.columns):
        df[msms_on] = df[msms_on].apply(ast.literal_eval)
    return PeakFrame(df) 
