## 通用版的Enrich

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection 


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



import requests
import pandas as pd

class RaMPAnalysis:
    def __init__(self, base_url="https://ramp.niaid.nih.gov/api/v1"):
        """
        初始化RaMP API分析工具
        
        :param base_url: RaMP API的基础URL
        """
        self.base_url = base_url
    
    def enrich_pathway(self, ids, id_type='kegg'):
        """
        进行代谢通路富集分析
        
        :param kegg_id: KEGG代谢物ID列表
        :return: 富集分析结果DataFrame
        """
        # 为每个KEGG ID添加 'kegg:' 前缀
        prefixed_kegg_ids = [f"{id_type}:{id}" for id in ids]
        
        # 准备请求载荷
        payload = {
            "analytes": prefixed_kegg_ids
        }
        
        # RaMP API的基础URL
        base_url = "https://rampdb.nih.gov/api/pathways-from-analytes"
        
        try:
            # 发送POST请求到富集分析端点
            response = requests.post(
                f"{base_url}/pathway/enrichment", 
                json=payload,
                headers={'Content-Type': 'application/json'}
            )
            
            # 检查响应
            response.raise_for_status()
            
            # 解析结果
            results = response.json()
            
            # 转换为DataFrame
            if results and 'pathwayEnrichment' in results:
                df = pd.DataFrame(results['pathwayEnrichment'])
                
                # 添加关键统计信息列
                df['adjusted_pvalue'] = -np.log10(df['pValue'])  # 转换p值
                df = df.sort_values('adjusted_pvalue', ascending=False)
                
                # 返回处理后的结果
                return df
            else:
                print("未找到富集分析结果")
                return pd.DataFrame()
        
        except requests.RequestException as e:
            print(f"API请求错误: {e}")
            return pd.DataFrame()
    
    def disease_enrichment_analysis(self, metabolite_ids=None, gene_ids=None):
        """
        进行疾病富集分析
        
        :param metabolite_ids: 代谢物ID列表
        :param gene_ids: 基因ID列表
        :return: 疾病分析结果DataFrame
        """
        # 准备请求载荷
        payload = {
            "metaboliteIds": metabolite_ids or [],
            "geneIds": gene_ids or []
        }
        
        # 发送POST请求
        try:
            response = requests.post(
                f"{self.base_url}/disease/enrichment", 
                json=payload,
                headers={'Content-Type': 'application/json'}
            )
            
            # 检查响应
            response.raise_for_status()
            
            # 解析结果
            results = response.json()
            
            # 转换为DataFrame
            if results and 'diseaseEnrichment' in results:
                df = pd.DataFrame(results['diseaseEnrichment'])
                
                # 添加关键统计信息列
                df['adjusted_pvalue'] = -np.log10(df['pValue'])  # 转换p值
                df = df.sort_values('adjusted_pvalue', ascending=False)
                
                return df
            else:
                print("未找到疾病富集分析结果")
                return pd.DataFrame()
        
        except requests.RequestException as e:
            print(f"API请求错误: {e}")
            return pd.DataFrame()
    
    def save_results(self, df, filename):
        """
        保存分析结果到CSV文件
        
        :param df: 分析结果DataFrame
        :param filename: 保存的文件名
        """
        df.to_csv(filename, index=False)
        print(f"结果已保存到 {filename}")

