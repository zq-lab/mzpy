'''
不是以mz对标为基础的富集算法，是代谢物名称匹配为基础的富集算法
'''



import certifi
import json
import os
import pandas as pd
import urllib3

class Enrichment(pd.DataFrame):
    # 继续继承该类，会导致子类的info函数不可用
    @property
    def _constructor(self):
        return self.__class__
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    @classmethod
    def load_kegg_feature(cls, org:str):
        fpath = os.path.dirname(os.path.abspath(__file__))
        data_fpath = os.path.join(fpath, "data", "kegg_features.pkl")
        df = pd.read_pickle(data_fpath)
        cdt1 = df['category'].isin(("enzyme", "module", "reaction"))
        cdt2 = df['kid'].str.startswith(org)
        df = df[cdt1 | cdt2]
        return cls(df)

    def get_features(self,feature:str):
        '''
        feature, column name of feature
        '''
        return set(self[feature].explode())
               
    def enrich(self, keys:set, feature:str, n_ft:str,
                func:str= 'fisher'):
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

    def validate(self, fdr = 0.01, nlimit = 20):
        '''
        满足fdr阈值的前nlimit条
        '''
        return self[self['_fdr_'] < fdr].head(nlimit)


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
            if not df.empty:
                df = df.sort_values(by='Pval_FDR', ascending =True)
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
        response = self.self.http.request('GET', url)

    def http_get(self, url):
        response = self.self.http.request('GET', url)
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