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
import re
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

from . import ms


### basic function for mzFrame
def concat(mzframe_list, msms_on='msms', ignore_index=False):
    for df in mzframe_list:
        df[msms_on] = df[msms_on].apply(ms.to_str)
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

        self.msms = ms.normalize(msms)   

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
        return ms.centroid(self.msms,
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
        return ms.match_mz(self.precursormz, mz, tol=tol, tol_rel=tol_rel, mode=mode)

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
            
 

class mzFrame(pd.DataFrame):
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
        self['msms'] = self['msms'].apply(lambda x: ms.centroid(x))


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
            self, mzFrame or pandas data frame object
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
                           ms_on='msms',
                           tol=(0.003, 0.005),
                           precursormz_compared=True,
                           similarity=0.99,
                           keep_first_on = None,
                           match_class='cpu'):
        '''
        drop duplicated msms

        param:
            mz_on, ms_on:  columns names of precursor mz and MSMS spectra
            tol:           tolerance for match
            smililarity:   similarity threshold for duplicates judgement. 
                            Based on the kernel density analysis of the metabolomics data from zebrafish,
                            the cutoff value for the similarity of MSMS matches should be set between 
                            0.85 and 0.92.
            keep_first_on: if None, retain the first one in the order of appearance.
                            if not None, sort by the column name specified in this parameter, 
                                and then keep the first one.
            match_class:   determine to use cpu or gpu edition for Match class

        return
            a data frame after deduplicates.
        '''
        if keep_first_on:
            df = self.sort_values(by=keep_first_on).copy().reset_index()
        else:
            df = self.copy().reset_index()

        if match_class == 'gpu':
            from .match_cp import MSList_cp
            msl = MSList_cp(df[mz_on], df[ms_on])
        else:
            from .match import MSList
            msl = MSList(df[mz_on], df[ms_on])
        
        score = msl.compute_similarity_self(tol, precursormz_compared)
        
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
                if ms.match_mz(mz, pcs.loc[j, 'mz']) == True:                    
                    df.loc[idx,'precursortype'] = pcs.loc[j, 'type']
                    break
        return df[df['precursortype'] != '']
    
    def match(self, que_mz, que_ms, mz_on='precursormz', ms_on='msms', match_class_type='cpu'):
        '''
        Calculate the similarity matrix between self ion ms (n) and que ms (m), 
            with the matrix having n rows and m columns.
        '''
        if match_class_type == 'gpu':
            from .match_cp import Match
        else:
            from .match import Match

        mat = Match(self[mz_on], self[ms_on], que_mz, que_ms)
        score_mx = mat.cosine_mx()
        return score_mx

    def fit(self,
            index_on,
            que_mz,
            que_ms,
            mz_on='precursormz',
            ms_on='msms',
            score_thd=0.85,
            stat_test='fisher',
            match_class_type='cpu'):
        '''
        score can be a vector or matrix

        param
            index_on, enrich targets
        '''
        if index_on not in self.columns:
            raise KeyError(f'{index_on}: unkown base column.')
        if len(que_mz) != len(que_ms):
            raise ValueError(f'The lengths of que_mz ({len(que_mz)}) and que_ms ({len(que_ms)})do not match.')
        
        score = self.match(que_mz, que_ms, mz_on=mz_on, ms_on=ms_on,
                           match_class_type=match_class_type)
        score = score.max(axis=1)
               
        hit_df = self[score > score_thd]
        if hit_df.shape[0] == 0:
            return
        
        ft = self[index_on].value_counts().to_frame(name='n_feature')
        hit = hit_df[index_on].value_counts().to_frame(name='n_hit')
        result_df = ft.join(hit, how='right')
        result_df['hit_%'] = 100 * result_df['n_hit'] / result_df['n_feature']
        result_df['hit_%'] = result_df['hit_%'].round(2)

        n_total = self.shape[0]
        n_total_hit = hit_df.shape[0]
        num_test_ft = len(que_mz)
        if stat_test == 'fisher':
            '''
            fisher table 2*2
                当前匹配数,   当前不匹配数
                其它匹配数，  其它不匹配数
            ref: https://www.statology.org/fishers-exact-test/
            '''
            _test = result_df[['n_feature', 'n_hit']].apply(lambda row: 
                                stats.fisher_exact([[row['n_hit'],
                                                    row['n_feature']-row['n_hit']],
                                                    [n_total_hit - row['n_hit'], 
                                                    n_total - row['n_feature']]]),
                            axis=1,
                            result_type='expand')
            pval = _test[1].tolist() 

        elif stat_test == 'hypergeom':
            if num_test_ft is None:
                raise ValueError(f'num_test_ft is required: {num_test_ft}')
            pval = result_df[['n_hit', 'n_feature']].apply(lambda row: 
                                                stats.hypergeom.sf(row['n_hit']-1, 
                                                                    n_total,
                                                                    row['n_feature'],
                                                                    num_test_ft), ## 这个二联表顺序有问题
                                                                    # 导致pcal全部为1
                axis=1)
        else:
            raise ValueError('inproper function No (func_no)')
        
        result_df['pval'] = pval
        _, fdr = fdrcorrection(pval)
        result_df['fdr'] = fdr
        result_df['-lgFDR'] = -np.log10(fdr)
        return result_df.sort_values(by='-lgFDR', ascending=False)
        
    def match_precursor_mz(self, mz, mz_on='precursormz', tol=0.003, tol_rel=5, mode='abs'):
        '''
        return precursor mz matched result
        '''
        condition = self[mz_on].apply(lambda x: ms.match_mz(x, mz, tol=tol, tol_rel=tol_rel, mode=mode))
        return self[condition]
    
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
                    if ms.match_mz(df.loc[idx, precursormz_on], df.loc[idx, msms_on][0][0],
                               tol=tol,
                               tol_rel=tol_rel,
                               mode=mode):
                        idx_to_drop.append(idx)                
            df = df.drop(idx_to_drop)
        
        df[msms_on] = df[msms_on].apply(lambda msms: ms.normalize(msms))
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
        peaks = ms.centroid(self[[mz_on, intensity_on]].values)
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

    def to_pickle(self, path, msms_on='msms', to_msms_str=True, *args, **kwargs):
        df = self.copy()
        if to_msms_str and (msms_on in df.columns):
            df[msms_on] = df[msms_on].apply(ms.to_str)            
        return super().to_pickle(path, *args, **kwargs)
    
    def to_sqlite3(self, tbl_name, conn, if_exists='replace', index=False, msms_on='msms'):
        df = self.copy()
        if msms_on in df.columns:
            df[msms_on] = df[msms_on].apply(ms.to_str) # 仅当msms是np数组时才有效
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
    mgf = mzFrame(data)
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
    df = mzFrame. read_msp(fpath, sep_ms2=sep_ms2)
    if extract_smiles:
        df['smiles'] = df['Comments'].str.extract('SMILES=(.*?)"')
    return df


def read_msd_ali(fname, washed=False):
    from .metab import Metab
    return Metab.read_msd_alignment(fname, washed=washed)


def read_msd_msp(mzFrame, fname, **kwargs):
    '''
    read peak list (msp format) exported from MS-Dial version 5.2 or higher
    peak height and peak area are transfered into relative value to the base peak.
    return:
        a mzFrame
    '''
    df = read_msp(fname, ** kwargs)
    df[['pkid', 'peak_height', 'peak_area']] = df['comment'].apply(_ex_intensity_).apply(pd.Series)
    basepk_heght = df['peak_height'].max()
    basepk_area  = df['peak_area'].max()
    df['peak_height'] = 100*df['peak_height']/basepk_heght
    df['peak_area']   = 100*df['peak_area']/basepk_area
    return df


def read_msp(fpath,
             sep_ms2='\t',
             rename: dict=None,
             comment=None,
             to_float: set = {'PRECURSORMZ','RETENTIONTIME', 'INTENSITY'},
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
    ions = [mzFrame._parse_msp_txt(block, sep_ms2=sep_ms2) for block in txt_blocks]
    df = mzFrame(ions)

    # 将指定列转换为浮点数类型
    for col in to_float.intersection(df.columns):
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df


def read_pickle(fname, msms_on='msms', force_msms=False):
    '''
    read pickle file of mzFrame
    param:
        fname, pickle file name
        msms_on, column name for MSMS
        force_msms, Whether to force the parsing conversion of MSMS, 
            pd.read_pickle seems to automatically convert list strings to arrays.

    由于模块名更改，在导入以前旧的模块名保存的pickle文件时，再度加载会找不到原模块名而报错。
    解决方案是加入模块别名：
    from mzpy import mzPandas as mpd
    sys.modules['mzpy.mzFrame'] = mpd
    '''
    df = pd.read_pickle(fname)
    if force_msms and (msms_on in df.columns):
        df[msms_on] = df[msms_on].apply(ast.literal_eval)
    return mzFrame(df)


def read_sql(sql, conn, msms_on='msms', force_msms=True):
    '''
    the method pd.read_sql does not automatically convert array strings to arrays by default, 
        so force_msms default is True (open).
    '''
    df = pd.read_sql(sql, conn)
    if force_msms and (msms_on in df.columns):
        df[msms_on] = df[msms_on].apply(ast.literal_eval)
    return mzFrame(df) 
