## 通用版的Enrich

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection 

# from plotnine import ggplot, aes, geom_point, scale_size, theme_bw, labs, xlab, ylab




class EnrichX(pd.DataFrame):
    '''
    Enrich基础表样式

    '''
    # 继续继承该类，会导致子类的info函数不可用
        
    @property
    def _constructor(self):
        return self.__class__    

    def check_columns(self):
        __colnames__ = ['id', 'name', 'features']
        if not all(col in self.columns for col in __colnames__):  
            raise ValueError(  
                f"DataFrame must contain columns {__colnames__}, "  
                f"but got {list(self.columns)}."  
            )   
     

    @classmethod
    def create(cls, df, target_id_on, target_name_on, feature_on):
        '''
        param:
            target_id_on, columns name of target id
            target_name_on, column name of target name
            features_on, column name of feature list 
        '''

        tmp = (  
            df.groupby([target_id_on, target_name_on])[feature_on]  
            .apply(set)  
            .reset_index(name="features")  
        ) 
        tmp.columns = ['id', 'name', 'features']
        tmp['n_ft'] = tmp['features'].apply(len)
        return cls(tmp) 

    @property
    def unique_features(self):
        return self['features'].explode().unique()  

    @property
    def nunique_features(self):
        return self['features'].explode().nunique()   


    def hypergeom(self, observed_features:list):
        '''
        param:
            target_id_on, columns name of target id
            target_name_on, column name of target name
            features_on, column name of feature list 
            observed_features, a list or set of features acquired from experiment or test 
        return:
            pd.DataFrame  
        '''
        self.check_columns()
        observed_features = set(observed_features) & set(self.unique_features)
        n_observed_features = len(observed_features)
        M = self.nunique_features # 特征总数
        enr = self.copy()
 
        enr['hits'] = enr['features'].apply(lambda x: set(x) & observed_features)
        enr['n_hits'] = enr['hits'].apply(len)
        enr =enr[enr['n_hits'] > 0]
        enr['n_expected'] = enr['n_ft'].apply(lambda x: x * n_observed_features / M)
            # 在随机抽取 n 个蛋白（不考虑任何生物学信息）时，预期平均“落在该通路”的蛋白数。
            # 若实测的命中数 远高于这个期望值，则暗示这条通路可能存在富集现象。
        enr['hit_percent'] = enr['n_hits'] / enr['n_ft']
       
        pval = enr.apply(lambda row: 
                            stats.hypergeom.sf(row['n_hits']-1, 
                                                M,
                                                row['n_ft'],
                                                n_observed_features),
                            axis=1) 
        enr['pval'] = pval  
        _, fdr, = fdrcorrection(pval)  
        enr['fdr'] = fdr
        enr['pFDR'] = -np.log10(fdr)
        return enr.sort_values(by='fdr', ascending=True)
    

    def fisher(self, observed_features:list):
        '''
        构建2x2列联表  
        |-------------|-----------|--------------|
        |             | 差异代谢物 | 非差异代谢物  |  
        |-------------|-----------|--------------|  
        | 属于该通路   | a         | c            |  
        | 不属于该通路 | b         | d            |  

        submit [[a, b], [c, d]] to Fisher function.
        '''
        self.check_columns()
        observed_features = set(observed_features) & set(self.unique_features)
        n_observed_features = len(observed_features)
        M = self.nunique_features # 特征总数
        nondem = set(self.unique_features) - observed_features      # 非差异特征（代谢物） 

        enr =self.copy()

        enr['hits'] = enr['features'].apply(lambda x: set(x) & observed_features)
        enr['n_hits'] = enr['hits'].apply(len) 
        enr = enr[enr['n_hits'] > 0]
        enr['n_expected'] = enr['n_ft'].apply(lambda x: x * n_observed_features / M)
            # 在随机抽取 n 个蛋白（不考虑任何生物学信息）时，预期平均“落在该通路”的蛋白数。
            # 若实测的命中数 远高于这个期望值，则暗示这条通路可能存在富集现象。
        enr['hit_percent'] = enr['n_hits'] / enr['n_ft']
        
        pval = enr.apply(lambda row: 
                stats.fisher_exact([[row['n_hits'],              # 差异代谢物中属于该通路的数量
                                    row['n_ft'] - row['n_hits']], # 差异代谢物中不属于该通路的数量 
                                    [len(nondem & set(row['features'])),  # 非差异代谢物中属于该通路的数量
                                    M - len(nondem & set(row['features']))] # 非差异代谢物中不属于该通路的数量
                                ])[1],
            axis=1)
    
        enr['pval'] = pval  
        _, fdr, = fdrcorrection(pval)  
        enr['fdr'] = fdr
        enr['pFDR'] = -np.log10(fdr)
        return enr.sort_values(by='fdr', ascending=True)


