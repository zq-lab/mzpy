import numpy as np
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection 

def enrich_df(df, n_total_feature,  n_total_hit, fdr=True, method='hypergeom'): 
    """  
    对特征进行Fisher精确检验  
    
    参数:  
        df: pandas DataFrame, 特征匹配结果  
            第1列: n_feature (特征总数)  
            第2列: n_hit (特征命中数)  
        n_total: 总样本数  
        n_total_hit: 总命中数  
        fdr: 是否进行假发现率(FDR)校正  
        method: Fisher检验的备择假设类型 ('two-sided', 'less', 'greater')  
    
    返回:  
        p值数组  
    """  
    if method not in ('fisher', 'hypergeom'):
        ValueError(f'Unknown test method {method}') 

    enr = df.copy()    
    if method == 'fisher':
        enr['pval'] = fisher(enr, n_total_feature, n_total_hit, fdr=fdr)
    else:
        enr['pval'] = hypergeom(enr, n_total_feature, n_total_hit, fdr=fdr)
    
    enr['ratio'] = enr['num_hit_features'] / enr['num_features']
    pval = enr['pval'].values.copy()
    pval[pval == 0] = np.min(pval[pval > 0]) / 10  # 0 替换为最小值正数的1/10
    enr['-log_p'] = -1 * np.log10(pval)
    enr['score'] = enr['ratio'] * enr['-log_p'] 
    enr = enr.sort_values(by='score', ascending=False)
    return enr



def fisher(df, n_total_feature,  n_total_hit, fdr=True, method='two-sided'):  
    """  
    对特征进行Fisher精确检验  
    
    参数:  
    df: pandas DataFrame, 特征匹配结果  
        第1列: n_feature (特征总数)  
        第2列: n_hit (特征命中数)  
    n_total: 总样本数  
    n_total_hit: 总命中数  
    fdr: 是否进行假发现率(FDR)校正  
    method: Fisher检验的备择假设类型 ('two-sided', 'less', 'greater')  
    
    返回:  
    p值数组  
    """  

    # 类型检查和转换  
    a = df.iloc[:, 1].values.astype(np.int64)  # 特征组命中数  
    b = (df.iloc[:, 0].values - a).astype(np.int64)  # 特征组未命中数  
    c = (n_total_hit - a).astype(np.int64)  # 总体命中数 - 特征组命中数  
    d = (n_total_feature - df.iloc[:, 0].values).astype(np.int64)  # 总特征数 - 每个特征组总数   

    if (b < 0).any():
        raise ValueError('negative values found in b.')
    if (c < 0).any():
        raise ValueError('negative values found in c.')
    if (d < 0).any():
        raise ValueError('negative values found in d.')

    # 使用 NumPy 的 vectorize 进行矢量化计算  
    vec_fisher = np.vectorize(  
        lambda a, b, c, d: stats.fisher_exact([[a, b], [c, d]], alternative=method)[1],  
        otypes=[float]  
    )  
    
    pval = vec_fisher(a, b, c, d)  

    # FDR校正  
    if fdr:  
        # 使用 Benjamini-Hochberg 方法  
        _, corrected_pval = fdrcorrection(pval, method='indep')  
        return corrected_pval  
    else:  
        return pval



def hypergeom(df, n_total_feature, n_total_hit, fdr=True, method='greater'):  
    """  
    对特征进行超几何分布检验  
    
    参数:  
    df: pandas DataFrame, 特征匹配结果  
        第1列: n_feature (特征总数)  
        第2列: n_hit (特征命中数)  
    n_total: 总样本数  
    n_total_hit: 总命中数  
    fdr: 是否进行假发现率(FDR)校正  
    method: 检验方法 ('greater', 'less', 'two-sided')  
    
    返回:  
    p值数组  
    """  
    # 类型检查和转换  
    n_feature = df.iloc[:, 0].values.astype(np.int64)  
    n_hit = df.iloc[:, 1].values.astype(np.int64)  

    # 使用 NumPy 的 vectorize 进行矢量化计算  
    def _calc_pval(n_feature, n_hit):  
        if method == 'greater':  
            # 大于等于观测命中数的概率  
            return stats.hypergeom.sf(n_hit - 1, n_total_feature, n_total_hit, n_feature)  
        elif method == 'less':  
            # 小于等于观测命中数的概率  
            return stats.hypergeom.cdf(n_hit, n_total_feature, n_total_hit, n_feature)  
        else:  # two-sided  
            # 计算双侧p值  
            # 找到对称的极端值  
            p_greater = stats.hypergeom.sf(n_hit - 1, n_total_feature, n_total_hit, n_feature)  
            p_less = stats.hypergeom.cdf(n_hit, n_total_feature, n_total_hit, n_feature)  
            return 2 * min(p_greater, p_less)  

    vec_hypergeom = np.vectorize(_calc_pval, otypes=[float])  
    pval = vec_hypergeom(n_feature, n_hit)  

    # FDR校正  
    if fdr:  
        # 使用 Benjamini-Hochberg 方法  
        _, corrected_pval = fdrcorrection(pval, method='indep')  
        return corrected_pval  
    else:  
        return pval 