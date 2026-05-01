import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection 

def enrich_df(df, n_total_feature,  n_total_hit, use_fdr=True, method='hypergeom'): 
    """  
    对特征进行Fisher精确检验  
    
    参数:  
        df: pandas DataFrame, 特征匹配结果  
            第1列: n_feature (特征总数)  
            第2列: n_hit (特征命中数)  
        n_total: 总样本数  
        n_total_hit: 总命中数  
        use_fdr: 是否进行假发现率(FDR)校正  
        method: Fisher检验的备择假设类型 ('two-sided', 'less', 'greater')  
    
    返回:  
        p值数组  
    """  
    if method not in ('fisher', 'hypergeom'):
        ValueError(f'Unknown test method {method}') 

    enr = df.copy()    
    if method == 'fisher':
        enr['pval'] = fisher(enr, n_total_feature, n_total_hit, use_fdr=use_fdr)
    else:
        enr['pval'] = hypergeom(enr, n_total_feature, n_total_hit, use_fdr=use_fdr)
    
    enr['ratio'] = enr.iloc[:, 1] / enr.iloc[:, 0]
    pval = enr['pval'].values.copy()
    pval[pval == 0] = np.min(pval[pval > 0]) / 10  # 0 替换为最小值正数的1/10
    enr['-log_p'] = -1 * np.log10(pval)
    enr['score'] = enr['ratio'] * enr['-log_p'] 
    enr = enr.sort_values(by='score', ascending=False)
    return enr



def fisher(df, n_total_feature,  n_total_hit, use_fdr=True, method='two-sided'):  
    """  
    对特征进行Fisher精确检验  
    
    参数:  
    df: pandas DataFrame, 特征匹配结果  
        第1列: n_feature (特征总数)  
        第2列: n_hit (特征命中数)  
    n_total: 总样本数  
    n_total_hit: 总命中数  
    use_fdr: 是否进行假发现率(FDR)校正  
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
    if use_fdr:  
        # 使用 Benjamini-Hochberg 方法  
        _, corrected_pval = fdrcorrection(pval, method='indep')  
        return corrected_pval  
    else:  
        return pval



def hypergeom(df, n_total_feature, n_total_hit, use_fdr=True, method='greater'):  
    """  
    对特征进行超几何分布检验  
    
    参数:  
    df: pandas DataFrame, 特征匹配结果  
        第1列: n_feature (特征总数)  
        第2列: n_hit (特征命中数)  
    n_total: 总样本数  
    n_total_hit: 总命中数  
    use_fdr: 是否进行假发现率(FDR)校正  
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
    if use_fdr:  
        # 使用 Benjamini-Hochberg 方法  
        _, corrected_pval = fdrcorrection(pval, method='indep')  
        return corrected_pval  
    else:  
        return pval 
    


def plsda(X, y, n_components=2):
    '''
    X, a matrix, rows are samples, columns are features (genes, proteins, or metabolites)   
    y, group factor (numpy 1d array)
    '''
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.preprocessing import OneHotEncoder
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_auc_score
  
    # 对 y 进行独热编码
    enc = OneHotEncoder(sparse_output=False, dtype=float)
    Y = enc.fit_transform(np.array(y).reshape(-1, 1))
    class_names = enc.categories_[0].tolist()

    # 拟合 PLS-DA 模型
    pls = PLSRegression(n_components=n_components, scale=False)
    pls.fit(X, Y)

    # 得分（样本在潜变量空间的坐标）
    T_scores = pls.x_scores_[:, :2]  # 选择前两个潜变量

    # 交叉验证
    # 计算每个类别的最小样本数
    min_samples_per_class = min(pd.Series(y).value_counts())
    n_splits = min(5, min_samples_per_class)  # 合理的折叠数
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    def fit_and_predict(Xtr, ytr, Xte):
        enc_cv = OneHotEncoder(sparse_output=False, dtype=float).fit(ytr.reshape(-1, 1))
        Ytr = enc_cv.transform(ytr.reshape(-1, 1))
        model = PLSRegression(n_components=n_components, scale=False).fit(Xtr, Ytr)
        
        Yte_hat = model.predict(Xte)
        cols = enc_cv.categories_[0].tolist()
        aligned = np.zeros((Yte_hat.shape[0], len(class_names)))
        for j, cname in enumerate(cols):
            aligned[:, class_names.index(cname)] = Yte_hat[:, j]
        return aligned

    # 生成折外预测
    Yhat_cv = np.zeros((X.shape[0], len(class_names)))
    for (tr, te) in cv.split(X, y):
        Yhat_cv[te] = fit_and_predict(X[tr], y[tr], X[te])

    y_cv = np.array(class_names)[np.argmax(Yhat_cv, axis=1)]
    
    # 计算评价指标
    acc = accuracy_score(y, y_cv)
    cm = confusion_matrix(y, y_cv, labels=class_names)
    try:
        auc_ovr = roc_auc_score(OneHotEncoder(sparse_output=False).fit_transform(np.array(y).reshape(-1, 1)),
                                Yhat_cv, average="macro", multi_class="ovr")
    except Exception:
        auc_ovr = np.nan

    # 打印结果
    print(f"CV accuracy (5-fold): {acc:.3f}, macro-AUROC (OvR): {auc_ovr:.3f}")
    print(pd.DataFrame(cm, index=class_names, columns=class_names))
    print(classification_report(y, y_cv, target_names=class_names))

    return pls, T_scores


def vip_scores_plsda(pls_model):
    """
    计算 PLS-DA 的 VIP 分数（Wold VIP）。
    - pls_model: 已拟合的 sklearn.cross_decomposition.PLSRegression
      要求：已拟合，n_components = 你使用的成分数
    返回：
      vip: (n_features,) 的一维数组
    """
    T = pls_model.x_scores_            # (n_samples, n_comp)
    W = pls_model.x_weights_           # (n_features, n_comp)
    Q = pls_model.y_loadings_          # (n_targets, n_comp)  注意：多分类时 n_targets > 1

    # 计算每个成分对 Y 的平方和贡献 SSY_k，形状 (n_comp,)
    # 对于多响应 Y，SSY_k = sum_i || t_k * q_{i,k} ||^2 = (t_k^T t_k) * sum_i q_{i,k}^2
    # 其中 t_k 是第 k 个 X-score（长度 n_samples），q_{i,k} 是第 k 个 Y-loading 的第 i 个分量
    tt = np.sum(T**2, axis=0)             # (n_comp,)
    qq = np.sum(Q**2, axis=0)             # (n_comp,)
    SSY = tt * qq                          # (n_comp,)

    # W 的每个成分的范数平方，用来归一化
    Wnorm2 = np.sum(W**2, axis=0)          # (n_comp,)

    p = W.shape[0]                         # n_features
    vip = np.zeros(p)
    for j in range(p):
        weights_j = (W[j, :]**2) / Wnorm2  # (n_comp,)
        vip[j] = np.sqrt(p * np.sum(SSY * weights_j) / np.sum(SSY))
    return vip
