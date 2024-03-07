'''
DEM analysis and enrichment
    DEM, Differential Expression Metabolites
    A special DataFrame with secondary column headings
    and its associated DEM analysis and drawing method

    For MultiIndex, slice(None) can be used as placeholder: df[(slice(None), 'aa'), :]

默认：Metabolite name的;分隔的第一个值是代谢物的关键索引值

'''
import certifi
import json
import numpy as np
import os
import pandas as pd
from scipy import stats
from sklearn.impute import KNNImputer
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.stats.multitest import fdrcorrection
import urllib3

from .plotfine import PlotFine

class Metab(pd.DataFrame): 
    _id_pattern = {'kegg': '(C\d{5})',
                   'hmdb': '(HMDB\d{7})'}
    plot = PlotFine()   # 这导致覆盖了父类的plot函数


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return self.__class__
    
    @classmethod
    def read_MSdial_alignment(cls, fpath:str, washed = False):
        '''
        fpath: file path of MSdial-exported txt file 
        '''
        df = pd.read_table(fpath, header = [0,4],
                        #    index_col = 0,
                            low_memory = False)
        if 'NA' in df.columns:
            del df['NA']
        cols = [['feature', it[1]]
                    if it[0].startswith('Unnamed') or it[0].startswith('Class')
                    else list(it)
                        for it in df.columns]
        df.columns = pd.MultiIndex.from_tuples(cols)
        df = cls(df)
        if washed:
            df = df.wash()
        return df
   
    def wash(self, group_by:str='kid', total_score=1.0):
        '''
        - Extract kid from names.    
        - Quantitative values of the same components will be combined and added
        param:
            total_score, cutoff value for total score
        '''
        df = self[self[('feature', 'MS/MS matched')]].copy()
        df = df.sort_values(by=('feature', 'Total score'), ascending=False)
        df = df.drop_duplicates(subset=[('feature', 'INCHIKEY')])
        df.index = df[('feature', 'Metabolite name')].str.split(';', expand=True)[0].to_list()
        # 暂时忽略零值问题
        return df
            
    def fill_missing(self, n_neighbors:int = 5):
        '''only for missing data, not zero data'''
        impKNN = KNNImputer(n_neighbors = n_neighbors)
        new_values = impKNN.fit_transform(self)        
        return self.__class__(new_values, columns=self.columns, index=self.index)
    
    def fill_zero(self):
        groups = self.groups
        min_values = self[groups].apply(lambda x: x[x != 0].min()).min()
        df = self.copy()
        df.loc[:, groups] = df[groups].replace(0, min_values/10)
        return df


    @property
    def groups(self):
        return self.columns.levels[0].drop('feature').to_list()

    def extract_id(self,
                   target = 'kegg',
                   col_name = ('feature', 'Metabolite name'),                   
                   as_index = False):
        if target not in self._id_pattern.keys():
            raise ValueError('target must be one of %s'%(self._id_pattern.keys()))
        ids = self[col_name].str.extract(self._id_pattern[target], expand=True)[0]
        ids.name = target + '_id'
        if as_index:
            self.index = ids
        else:
            return ids
    
    def pca(self, groups:list = None, labeled=False, palette='Set1', save_to:str = None):
        df = self.fill_zero()
        if groups is None:
            groups = self.groups
        data = np.log10(df[groups].T)
        if labeled:
            labels = data.index.get_level_values(level=1)
        else:
            labels = None
        data.index = data.index.get_level_values(level=0)
        groups = pd.Categorical(data.index)
        return self.plot.pca(data, groups=groups, labels = labels, palette=palette,save_to=save_to)

    def vs(self, scheme:str,
            fc:float=1.5, p:float=0.05,
            palette = 'Set1', save_to:str = None):
        '''calculate g1/g2
        ---------------------------------------
        Definition standard of differential metabolites：
            Metabolites, mappable to KEGG or HMDB IDs, 
            that had a fold-change greater than +/− 1.5 
            with an FDR adjusted p-value <0.05
            ref: MEtabolites, 2018, https://www.mdpi.com/2218-1989/8/1/16
        parameters:
        -----------------
            scheme:  calculation scheme, for example, G1/G2
            fc, threshold value of fold change
            p, threhold value of p-value
            save_to, where to save volcano plot figure            
        return:
        ----------------
            return None
            - the vocano plot will be saved if set save_to
            - calculation results svaed into self data frame with scheme 'g1/g2'
        '''
        # scheme = f'_{scheme}_'
        g1, g2 = scheme.split('/')
        if (g1 not in self.columns) or (g2 not in self.columns):
            raise ValueError('scheme is imporper. It contains incorrect group name.')
        if scheme not in self.columns:
            # to avoide repeaded calculation
            # calc log2FC and q values
            avg_g1 = self[g1].apply(np.mean, axis=1)
            avg_g2 = self[g2].apply(np.mean, axis=1)
            # 找到均值中0以外的最小值
            min1 = avg_g1[avg_g1 != 0].min()
            min2 = avg_g2[avg_g2 != 0].min()
            # 均值的最小值的1/10用于替换0值
            # Avoid dividing by 0 in the next step
            min_sub = min(min1/10, min2/10)
            avg_g1 = avg_g1.replace(0, min_sub) 
            avg_g2 = avg_g2.replace(0, min_sub)

            self[(scheme, 'log2FC')] = np.log2(avg_g1/avg_g2)
            _, pval = stats.ttest_ind(self[g1], self[g2], axis = 1)
                # g1或g2如果有恒定的值，test会给出一个runtime warning
                # g1和g2完全相等，则pval为nan
            pval = np.where(np.isnan(pval), min(pval)/10, pval)
            self[(scheme,'pval')] = pval
            self[(scheme, '_fdr_')] = multipletests(pval, method='fdr_bh')[1]
            self[(scheme, '_pFDR_')] = self[(scheme, '_fdr_')].apply(lambda x: -np.log10(x))
        for i in self.index:
            # append trend schemes: up, dn (down), no
            if self.loc[i, (scheme, '_pFDR_')] > -np.log10(p):
                if self.loc[i, (scheme, 'log2FC')] > np.log2(fc):
                    self.loc[i, (scheme, 'trend')] = 'up'
                elif self.loc[i, (scheme, 'log2FC')] < -np.log2(fc):
                    self.loc[i, (scheme, 'trend')] = 'dn'
                else:
                    self.loc[i, (scheme, 'trend')] = 'no'
            else:
                self.loc[i, (scheme, 'trend')] = 'no'
        # ploting vocano digram
        plot = self.plot.volcano(self[scheme], x='log2FC', y='_pFDR_', fill='trend',
                                    xcut = np.log2(fc), ycut = -np.log10(p),
                                    title = scheme,
                                    palette = palette,
                                    save_to=save_to)
        return plot
        
    
    def vs3(self, groups:list, pattern:str,
            fc:float=1.5, p:float=0.05,
            key:str = None,
            palette = 'Set1',
            save_to = None):
        '''
        obtaine differential expressed metabolites (dem) form the existing vs groups (vs1 and vs2) according to
            the specified pattern, and plot venn diagram.

        param:
            vs1 and vs2, vs result yield from vs function
            pattern, anti, syn or var, Trends pattern among the three groups. 
                anti, with opposite trends;
                syn, with the same trends,;
                var, with the same or opposite trends
            key, if it is None, the row index will be used as identifiers. Or the specific column will be used as identifiers.
            plotted,  if it is true, a associated Venn diagram will be printed.
            save_to, save Venn diagrame to a specific path. If plotted is False, the item will be ignored
        '''
        if not isinstance(groups, list):
            raise TypeError('groups must be a list with 3 elements.')
        if len(groups) != 3:
            raise TypeError('groups must contains 3 elemnets.')
        if not set(groups).issubset(self.groups):
            raise ValueError('Unknown name in groups.')
        
        if save_to is not None:
            # label = '_'.join(groups)
            # save_to = save_to + label +'/'
            if not os.path.exists(save_to):
                os.makedirs(save_to)
            f = save_to + '%d_%s.svg'
        vs1 = f'{groups[1]}/{groups[0]}'
        vs2 = f'{groups[2]}/{groups[1]}'

        self.pca(groups, save_to= f%(1, 'pca'))

        self.vs(vs1, fc = fc, p = p, palette=palette,
                save_to=f%(2, vs1.replace('/', ''))
            )
        self.vs(vs2, fc = fc, p = p, palette=palette,
                save_to=f%(3, vs2.replace('/', ''))
            )
        
        if key:
            a = set(self.loc[self[(vs1, 'trend')] == 'up', ('feature', key)].values)
            b = set(self.loc[self[(vs1, 'trend')] == 'dn', ('feature', key)].values)
            c = set(self.loc[self[(vs2, 'trend')] == 'up', ('feature', key)].values)
            d = set(self.loc[self[(vs2, 'trend')] == 'dn', ('feature', key)].values)
        else:
            a = set(self[self[(vs1, 'trend')] == 'up'].index.values)
            b = set(self[self[(vs1, 'trend')] == 'dn'].index.values)
            c = set(self[self[(vs2, 'trend')] == 'up'].index.values)
            d = set(self[self[(vs2, 'trend')] == 'dn'].index.values)
        # gather data for venn plot
        if pattern in ('anti', 'syn'):
            data = {vs1+' up'  : a,
                    vs1+' down': b,
                    vs2+' up'  : c,
                    vs2+' down': d}                        
        elif pattern == 'var':
            data = {vs1: set(a|b),
                    vs2: set(c|d)}
        else:
            raise ValueError('The pattern should be one of  ("anti", "syn", "var")')
        if save_to:
            self.plot.venn(data, save_to = f%(3, 'Venn'))
        else:
            self.plot.venn(data)

        # gather dem
        if pattern == 'anti':
            dem = a.intersection(d) | b.intersection(c)
        elif pattern == 'syn':
            dem = a.intersection(c) | b.intersection(d)
        elif pattern == 'var':
            dem = (a|b).intersection(c|d)
        
        with open(save_to + '4_dem.txt', 'w') as txt_file:
            ss = '\n'.join(dem)
            txt_file.write(ss)
        
        self.to_csv(save_to + '/0_full_data.tsv', sep = '\t', index=False)
        return dem
        


class Enrichment(pd.DataFrame):
    # 继续继承该类，会导致子类的info函数不可用
    @property
    def _constructor(self):
        return self.__class__
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    

    def get_features(self,feature:str):
        '''
        feature, column name of feature
        '''
        return set(self[feature].explode())
               
    def enrich(self, keys:set, feature:str, n_ft:str,
                func:str = 'fisher'):
        '''
        tutorial for fisher test:
            https://www.statology.org/fishers-exact-test/
            https://github.com/zqfang/GSEApy/blob/master/gseapy/stats.py

        param:
            feature, column name of feature list (compounds)
            n_ft,    column name of feature number (n_compounds)
            func, 'fisher' or 'hypergeom'
        Note:
            For disease, Fisher test are more suitable
            For metabolism pathway, hypergeom may be more suitable.
        '''
        print('preparing data set')
        keys = set(keys)
        n_keys = len(keys)
        df = self[self[n_ft] > 0].copy()
        df['_match_'] = df[feature].apply(lambda x: set(x).intersection(keys))
        df['_n_match_'] = df['_match_'].apply(len)
        df['impact'] = df['_n_match_']/df[n_ft]
        n_total_features = len(df.get_features(feature)) 
        n_total_match = len(df.get_features(feature).intersection(keys))
        df = df[df['_n_match_'] > 0]
        print('enriching...')
        if func == 'fisher':
            '''
            fisher table 2*2
                当前匹配数,   当前不匹配数
                其它匹配数，  其它不匹配数
            ref: https://www.statology.org/fishers-exact-test/
            '''
            df[['_odds_ratio_', '_pval_']] = df[[n_ft, '_n_match_']].apply(lambda row: 
                stats.fisher_exact([[row['_n_match_'], row[n_ft]-row['_n_match_']],
                                    [n_total_match-row['_n_match_'], 
                                     n_total_features-row[n_ft]]]),
                axis=1,
                result_type='expand')
            # df['_odds_ratio_'] = odds_ratio
        elif func == 'hypergeom':
            '''
            ref: 
            '''
            df['_pval_'] = df[['_n_match_', n_ft]].apply(lambda row: 
                stats.hypergeom.sf(row['_n_match_'], 
                                    n_total_features,
                                    row[n_ft],
                                    n_keys),
                axis=1)
        else:
            raise ValueError('inproper function No (func_no)')       
        # df['_fdr_'] = multipletests(df['_pval_'].values, method='fdr_bh')[1]
        # 该方法与下面的函数计算结果一致，但multipletests的可选择方法更多
        _, df['_fdr_'], = fdrcorrection(df['_pval_'].values)
        df['_pFDR_'] = -np.log10(df['_fdr_'])
        df.sort_values(by='_fdr_', ascending=True, inplace=True) 
        return df


class RaMP():
    '''It can only be enriched to the pathway, not to the smaller nodes'''
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
                                ca_certs=certifi.where())
    url_root = url = 'https://rampdb.nih.gov/api/'                      

    def fisher_test(self, dem:list, id_source:str = 'kegg'):
        '''
        dem The dem id style needs to meet the requirements of RaMP
                https://rampdb.nih.gov/api, such like:
                    hmdb:HMDB0000201, chemspider:10026, CAS:5657-19-2, kegg:C00078
        
        param:
            id_source, ('hmdb','kegg', 'pubchem')
        '''
        dem = [f'{id_source.lower()}:{id}' for id in dem]
        url = self.url_root + 'combined-fisher-test'
        encoded_data = json.dumps({"analytes":  dem }).encode('utf-8')
        info = 'Please wait for the connection and acquire data from to RaMP official website'
        info += '\n due to the current network status.'
        print(info)
        response = self.http.request('POST', url,
                            body = encoded_data,  
                            headers = {'accept': 'application/json',
                                       'Content-Type': 'application/json'})
        if response.status == 200:
            result = response.data.decode('utf-8')
            rdata = json.loads(result)
            df = pd.DataFrame(rdata['data']['fishresults'])
            # possible empty df
            # df = df.sort_values(by='Pval_FDR', ascending =True)
            return df
        else:
            print('Response Error: %d'%response.status)

    @property
    def version(self):
        url = self.url_root + 'source_versions'
        response = self.http.request('GET', url)
        return self._post_process(response)
    
    @property
    def entity_counts(self):
        url = self.url_root + 'entity_counts'
        response = self.http.request('GET', url)

    def http_get(self, url):
        response = self.http.request('GET', url)
        return self._post_process(response)
    
    def http_post(self, url, ids:list):
        encoded_data = json.dumps({"analytes":  ids }).encode('utf-8')
        response = self.http.request('POST', url,
                            body = encoded_data,  
                            headers = {
                                'accept': 'application/json',
                                'Content-Type': 'application/json'})
        return self._post_process(response)

    def _post_process(self, response:urllib3.response.HTTPResponse)->pd.DataFrame:
        # return 
        if response.status == 200:
            result = response.data.decode('utf-8')
            rdata = json.loads(result)
            return  rdata['data']
        else:
            print('Response Error: %d'%response.status)