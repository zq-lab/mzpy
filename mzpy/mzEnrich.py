import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection


## 是否要继承mzFrame类？
class Enrich(pd.DataFrame):
    '''
    enrich base table structure:
        label 1, mz11, ms11
        label 1, mz12, ms12
        ...,     ..., ...
        label 1, ..., ...
        label 2, mz21, ms21
        ...
    '''
     
    @property
    def _constructor(self):
        return self.__class__
        # 使用pd.concat之后仍能保持子类类型
        # 使用pd.merge之后仍能保持子类类型
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def match(self, que_mz, que_ms, mz_on='precursormz', ms_on='msms', match_class_type='cpu'):
        if match_class_type == 'gpu':
            from .mzMatch_cp import Match
        else:
            from .mzMatch import Match

        mat = Match(self.base[mz_on], self.base[ms_on], que_mz, que_ms)
        score_mx = mat.cosine_mx()
        return score_mx
    
    def fit(self,
            index_on,
            que_mz,
            que_ms,
            mz_on='precursormz',
            ms_on='msms',
            score_thd=0.85,
            stat_test='fisher'):
        '''
        score can be a vector or matrix

        param
            index_on, enrich targets
        '''
        if index_on not in self.base.columns:
            raise KeyError(f'{index_on}: unkown base column.')
        if len(que_mz) != len(que_ms):
            raise ValueError(f'The lengths of que_mz ({len(que_mz)}) and que_ms ({len(que_ms)})do not match.')
        
        score = self.match(que_mz, que_ms, mz_on=mz_on, ms_on=ms_on)
        score = score.max(axis=1)
               
        hit_df = self.base[score > score_thd]
        if hit_df.shape[0] == 0:
            return
        
        ft = self.base[index_on].value_counts().to_frame(name='n_feature')
        hit = hit_df[index_on].value_counts().to_frame(name='n_hit')
        result_df = ft.join(hit, how='right')
        result_df['hit_%'] = 100 * result_df['n_hit'] / result_df['n_feature']
        result_df['hit_%'] = result_df['hit_%'].round(2)

        n_total = self.base.shape[0]
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
                                                                    num_test_ft),
                axis=1)
        else:
            raise ValueError('inproper function No (func_no)')
        
        result_df['pval'] = pval
        _, fdr = fdrcorrection(pval)
        result_df['fdr'] = fdr
        result_df['-lgFDR'] = -np.log10(fdr)
        return result_df.sort_values(by='-lgFDR', ascending=False)
    
    