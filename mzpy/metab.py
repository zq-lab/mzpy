'''
DEM analysis and enrichment
    DEM, Differential Expression Metabolites
    A special DataFrame with secondary column headings
    and its associated DEM analysis and drawing method

    For MultiIndex, slice(None) can be used as placeholder: df[(slice(None), 'aa'), :]

默认：Metabolite name的;分隔的第一个值是代谢物的关键索引值

'''
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import KNNImputer
from statsmodels.sandbox.stats.multicomp import multipletests

from .mzPandas import mzFrame
from .mzPlot import Plot
from . import ms

_id_pattern_ = {'kegg': r'(C\d{5})',
                'hmdb': r'(HMDB\d{7})'}

mzplt = Plot()

class Metab(pd.DataFrame):  

    @property
    def _constructor(self):
        return self.__class__
    
    def drop_duplicated_ms(self, 
                           mz_on='Average Mz',
                           ms_on='msms',
                           tol=(0.003, 0.005),
                           similarity=0.99,
                           keep_first_on = 'S/N average',
                           match_class='cpu'):
        '''
        drop duplicated msms or metabolites identified
        Note: 这里没有按MS/MS matched和INCHIKEY去重复。鉴定物的重复性是比较复杂的问题：
            首先，鉴定过程中，在比对MS相似度时，有阈值设定的问题，鉴定为同一个物质的离子，有可能并非时相同的离子，因此不能直接消除
            再者，列表中给INCHIKEY行，有的仅仅时一级离子匹配
            总之，要如何去重复，要根据实际的情况处理


        param:
            mz_on, ms_on:  columns names of precursor mz and MSMS spectra
            tol:           tolerance for match
            smililarity:   similarity threshold for duplicates judgement
            keep_first_on: if None, retain the first one in the order of appearance.
                            if not None, sort by the column name specified in this parameter, 
                                and then keep the first one.
            match_class:   determine to use cpu or gpu edition for Match class

        return
            a data frame after deduplicates.
        '''
        # 按转为mzFrame后按msms去重复，获得去重复后的行索引
        mdf = mzFrame(self['_'])
        mdf = mdf.sort_values(by=keep_first_on, ascending=False)        
        mdf = mdf.drop_duplicated_ms(mz_on=mz_on,
                                     ms_on=ms_on,
                                     tol=tol,
                                     similarity=similarity,
                                     keep_first_on=keep_first_on,
                                     match_class=match_class)
        
        return self.loc[mdf.index]

    def extract_id(self,
                   target = 'kegg',
                   col_name = ('_', 'Metabolite name'),                   
                   as_index = False):
        if target not in _id_pattern_.keys():
            raise ValueError('target must be one of %s'%(_id_pattern_.keys()))
        ids = self[col_name].str.extract(_id_pattern_[target], expand=True)[0]
        ids.name = target + '_id'
        if as_index:
            self.index = ids
        else:
            return ids
               
    def fill_missing(self, n_neighbors:int = 5):
        '''only for missing data, not zero data'''
        impKNN = KNNImputer(n_neighbors = n_neighbors)
        new_values = impKNN.fit_transform(self)        
        return self.__class__(new_values, columns=self.columns, index=self.index)
    
    def fill_quantum_zero(self, replace_value=None):
        '''
        清洗定量数据
        '''
        df = self.copy()
        groups = df.groups

        if replace_value is None:
            qa = df.quantum
            replace_value = qa[qa > 0].min().min() / 10

        df.loc[:, groups] = df.loc[:, groups].replace(0, replace_value)
        #确保定量子表都是numeric类型
        # df.loc[:, groups] = df.loc[:, groups].apply(pd.to_numeric, errors='coerce')

        return df
    
    def get_factor(self, groups):
        cols_0 = self.columns.get_level_values(0)
        out_idx = set(groups) - set(cols_0)
        if len(out_idx) > 0:
            raise KeyError(f'not found keys: {out_idx}')
        else:
            return self[groups].columns.get_level_values(level=0)

    @property
    def groups(self):
        return self.columns.levels[0].drop('_').to_list()
    
    def hstack(self, _cols, groups, drop_col_levels=False):
        '''
        按选定的代谢物信息列和分组数据列选取子表后横向拼接        
        '''
        df_a = self.loc[:,('_', _cols)]
        df_b = self[groups]
        df = pd.concat([df_a, df_b], axis=1)
        if drop_col_levels:
            df.columns = df.columns.droplevel(0) 
            df = pd.DataFrame(df) # 没有了上级信息，这是要转换为普通的数据框

        return df 
    
    def log2FC(self, nume, deno):
        '''
        calculate log2(fold change): numer / deno
        param:
            numer, numerator group name
            deno, denominatorgroup name
        return:
            a vector
        '''

        avg_nume = self[nume].apply(np.mean, axis=1)
        avg_deno = self[deno].apply(np.mean, axis=1)

        # 找到均值中0以外的最小值,
        min_nume = avg_nume[avg_nume != 0].min()
        min_deno = avg_deno[avg_deno != 0].min()
        # Avoid dividing by 0 in the next step
        avg_nume = avg_nume.replace(0, min_nume) 
        avg_deno = avg_deno.replace(0, min_deno)

        return np.log2(avg_nume/avg_deno)        


    def pca(self, groups: list = None, labeled=False, palette='Set1', save_to: str = None):       
        # 如果未指定 groups，使用默认的 groups
        if groups is None:
            groups = self.groups

        # data = self.wash_quantum()
        data = self[groups].T.copy()
        
        # 对数据取对数
        data = np.log10(data)
        
        # 处理标签
        if labeled:
            labels = data.index.get_level_values(level=1)
        else:
            labels = None
        
        # 设置索引
        data.index = data.index.get_level_values(level=0)
        groups = pd.Categorical(data.index)
        
        # 调用绘图函数
        return mzplt.pca(data, groups=groups, labels=labels, palette=palette, save_to=save_to)
    
    @property
    def quantum(self, groups=None):
        '''
        return the quantity data frame
        '''
        if not groups:
            groups = self.groups
        return self[groups]
    
    def quantum_melt(self, id_vars, groups, value_name='peak area'):
        '''
        定量表改为长表，转后成长表后根据id注释代谢物的名称会更方便
        param
            id_vars, id for metabolites
            groups
        retunr:
            long table contains three columns: id, group, peak area (value_name)
        '''
        df = self.wash_metabolites().copy()
        df = self.hstack(_cols = id_vars, groups= groups, drop_col_levels=True)
        factors = self.get_factor(groups)
        df.set_index(id_vars, inplace=True)
        df.columns = factors
        df = df.reset_index()
        data = df.melt(id_vars=id_vars, value_name=value_name)
        data.columns = ['id', 'group', 'peak area']
        return data

    @property
    def scores(self):
        '''
        describe of all scores if MS matching.
        '''
        return self['_'][['m/z similarity',
           'Simple dot product', 'Weighted dot product', 'Reverse dot product',
           'Matched peaks count', 'Matched peaks percentage', 'Total score']].describe()


    def spearman(self, groups=None, corr_cut=0.7, p_cut=0.05, fdr_corr=True, save_to=None):
        '''
        param:
            groups, list, groups for Spearman's test
            corr, cutoff for correlation. 
                When |corr| > 0.7, it is typically considered to be strongly correlated 
                and suitable as a high-confidence marker.
            p, cutoff for p values
            fdr_corr whether to be corrected by FDR
        return:
            a data frame containing Spearman's results derived form self
        '''
        from scipy.stats import spearmanr

        if groups is None:
            groups = df.groups
            df = self.copy()
        else:
            df = self[['_'] + groups].copy()

        corr = []
        pval = []
        factor = df.get_factor(groups)
        # need check const values
        for i in df.index:
            correlation, p = spearmanr(df.loc[i, groups].values, factor)
            corr.append(float(correlation))
            pval.append(float(p))

        correlation = np.asarray(correlation)
        correlation[np.isnan(correlation)] = 0
        pval = np.asarray(pval)
        pval[np.isnan(pval)] = 1.0   
             
        if fdr_corr:
            pval = multipletests(pval, method='fdr_bh')[1]
        
        df[('spearman', 'corr')] = corr
        df[('spearman', 'pval')] = pval
        df[('spearman', '-lg_pval')] = -np.log10(pval)

        for i in df.index:
            if (df.loc[i, ('spearman', 'corr')] >= corr_cut) and \
            (df.loc[i, ('spearman', 'pval')] < p_cut):
                df.loc[i, ('spearman', 'monot')] = 'up'
            elif (df.loc[i, ('spearman', 'corr')] <= -1 * corr_cut) and \
                 (df.loc[i, ('spearman', 'pval')] < p_cut):
                df.loc[i, ('spearman', 'monot')] = 'dn'
            else:
                df.loc[i, ('spearman', 'monot')] = 'no'

        return df, mzplt.volcano(df['spearman'],
                                 x = 'corr',
                                 y = '-lg_pval',
                                 fill = 'monot',
                                 title = 'Spearman test',
                                 xlab = r'$r_s$',
                                 ylab = r'-$\log_{10}(\mathrm{p\text{-}value})$',
                                 xcut = corr_cut,
                                 ycut = -np.log10(p_cut),
                                 save_to=save_to)


    def trio(self, vs1, vs2,
            pattern:str,
            fc:float=1.5, p:float=0.05,
            key:str = None,
            palette = 'Set1',
            save_to = None):
        '''
        obtaine differential expressed metabolites (dem) form the existing vs groups (vs1 and vs2) according to
            the specified pattern, and plot venn diagram.
            vs must be run befor this fuction

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

        groups = set(vs1.split('/') + vs2.split('/'))
        if len(groups) != 3:
            raise TypeError('Only 3 groups can be accepted.')
        if not groups <= set(self.groups):
            raise ValueError('Unknown group name. Please Check groups names in params vs1 or vs2')        
      
        if key:
            a = set(self.loc[self[(vs1, 'regulation')] == 'up', ('_', key)].values)
            b = set(self.loc[self[(vs1, 'regulation')] == 'dn', ('_', key)].values)
            c = set(self.loc[self[(vs2, 'regulation')] == 'up', ('_', key)].values)
            d = set(self.loc[self[(vs2, 'regulation')] == 'dn', ('_', key)].values)
        else:
            a = set(self[self[(vs1, 'regulation')] == 'up'].index.values)
            b = set(self[self[(vs1, 'regulation')] == 'dn'].index.values)
            c = set(self[self[(vs2, 'regulation')] == 'up'].index.values)
            d = set(self[self[(vs2, 'regulation')] == 'dn'].index.values)
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
            mzplt.venn(data, save_to = save_to)
        else:
            mzplt.venn(data)

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
        
        self.to_csv(save_to + '0_full_data.tsv', sep = '\t', index=False)
        return dem

    def plsda(self, groups:list=None,palette='Set1', save_to:str=None):
        # 准备数据
        from sklearn.cross_decomposition import PLSRegression
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        if groups is None:
            groups = self.groups
        ft = self[groups]
        X = ft.T.values
        y = ft.columns.get_level_values(level=0).values
        y_std, labels = pd.factorize(y)

        # 数据预处理
        scaler = StandardScaler()
        X_std = scaler.fit_transform(X)

        pls_da = PLSRegression(n_components=2)
        X_plsda = pls_da.fit_transform(X_std, y_std)[0]
        lr = LogisticRegression(multi_class='ovr', random_state=1, solver='lbfgs')
        lr = lr.fit(X_plsda, y_std)
        # 绘制决策区域
        mzplt.decision_regions(X_plsda, y_std, labels, classifier=lr, cmap=palette,save_to=save_to)

        # 获取模型的加载向量
        loading_vectors = pls_da.x_weights_

        # 计算每个变量在每个主成分上的贡献度
        contributions = np.sum(loading_vectors**2, axis=1)

        # 计算 VIP 值
        t = pls_da.x_scores_
        w = pls_da.x_weights_
        q = pls_da.y_loadings_

        p, h = w.shape
        vips = np.zeros((p,))

        s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
        total_s = np.sum(s)

        for i in range(p):
            weight = np.array([w[i, j] ** 2 * s[j] for j in range(h)])
            vips[i] = np.sqrt(p * np.sum(weight) / total_s)
        return vips
    

        
    def ttest(self, g1, g2, fdr_corr=True):
        _, pval = stats.ttest_ind(self[g1].values, self[g2].values, axis = 1) 
        # g1或g2如果有恒定的值，test会给出一个runtime warning
        
        # g1和g2完全相等，则pval为nan
        pval = np.where(np.isnan(pval), np.nanmin(pval) / 10, pval)
        indeces = np.where(np.isnan(pval))
        if len(indeces) > 0 :
            print('identical or constant values are detected in lines:')
            print(self.index[indeces])
            print('These empty p-values will be replaced with 1/10 minimum p-value, and the calculation will continue.')
        
        pval = np.asarray(pval)
        pval[np.isnan(pval)] = 1.0

        if fdr_corr:
            return multipletests(pval, method='fdr_bh')[1]
        else:
            return pval

    def vs(self,
           nume,
           deno,
           fc:float=1.5,
           p:float=0.05,
           basic_info_on=['Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
                          'Adduct type',  'Formula',         'Ontology',    'INCHIKEY',
                           'SMILES',      'Total score'],
           ms_on = 'msms',
           palette = 'Set1',
           save_fig_to=None):
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
            vip_on: column name of PLS-DA vip values.
                If it is defined, p will be ignored.
            save_to, where to save volcano plot figure            
        return:
        ----------------
            return None
            - the vocano plot will be saved if set save_to
            - calculation results svaed into self data frame with scheme 'g1/g2'
        '''
        log2FC = self.log2FC(nume=nume, deno=deno)
        fdr_p = self.ttest(nume, deno)


        if basic_info_on:
            if ms_on and (ms_on in self['_'].columns):
                # 自动续加ms列
                basic_info_on.append(ms_on)
            df = self['_'][list(set(basic_info_on))].copy()
        else:
            df = pd.DataFrame()        


        df['log2FC'] = log2FC
        df['fdr_p'] = fdr_p
        df['-log10 (FDR p-val)'] = -1 * np.log10(fdr_p)
        df['regulation'] = df['log2FC'].apply(
                lambda x: 'up' if x > np.log2(fc) else (
                            'dn' if x < -np.log2(fc) else 'no')
        )
        ## 其次，所有p值不达标的均改为no
        df['regulation'] = df[['regulation', 'fdr_p']].apply(lambda row: 'no' if row['fdr_p'] >= p \
                                                   else row['regulation'],
                                                   axis=1)

        # ploting vocano digram
        plot = mzplt.volcano(df, x='log2FC', y='-log10 (FDR p-val)', fill='regulation',
                               xcut = np.log2(fc),
                               ycut = -np.log10(p),
                               title = f'{nume} / {deno}',
                               palette = palette,
                               save_to=save_fig_to)
        return df, plot        

    def wash_metabolites(self, total_score=1.0, keep_first_by='INCHIKEY'):
        '''  
        drop off unknown ions ans drop duplicated metabolites.
        param:
            total_score, cutoff value for total score
        '''
        df = self[self[('_', 'MS/MS matched')]].copy()
        df = df[df[('_', 'Total score')] > total_score]
        df = df.sort_values(by=('_', keep_first_by), ascending=False)
        df = df.drop_duplicates(subset=[('_', 'INCHIKEY')])
        df.index = df[('_', 'Metabolite name')].str.split(';', expand=True)[0].to_list()
        # 没有处理零值问题
        return df
    


    

def read_msd_ali(fpath:str, drop_null_ms=True):
    '''
    fpath: file path of MSdial-exported txt file 
    '''
    df = pd.read_table(fpath, header = [0,4],
                    #    index_col = 0,
                        low_memory = False)
    if 'NA' in df.columns:
        del df['NA']
    cols = [['_', it[1]]
                if it[0].startswith('Unnamed') or it[0].startswith('Class')
                else list(it)
                    for it in df.columns]
    df.columns = pd.MultiIndex.from_tuples(cols)
    df = Metab(df)

    if drop_null_ms:
        df = df[df[('_', 'MS/MS spectrum')].notnull()].reset_index(drop=True)
        df[('_', 'msms')] = df[('_', 'MS/MS spectrum')].apply(lambda x: ms.to_array(x, sep1=' ', sep2=':'))
        df.drop(columns=('_', 'MS/MS spectrum'))

    df = df.set_index(('_', 'Alignment ID'))
    # 重命名行索引为 'Alignment ID'
    df = df.rename_axis('Alignment ID')

    return df


