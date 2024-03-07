'''
HMDB xml data sheet processor
'''
import pandas as pd
from tqdm import tqdm
import xml.etree.ElementTree as ET
from .metab import Enrichment

class HmdbElement():
    '''
    采用了封装而非继承
    继承的问题太多
    '''
    tag_root = '{http://www.hmdb.ca}'

    def __init__(self, element:ET.Element):
        self._element = element

    @classmethod
    def load(cls, fpath):
        tree = ET.parse(fpath)
        root = tree.getroot()
        return cls(root)

    def __len__(self):
        return self._element.__len__()

    def __getitem__(self, index):
        return self.__class__(self._element[index])

    def __repr__(self):
        return f'<HmdbElement {self._element.tag}>'

    @property
    def tag(self):
        return self._element.tag
    
    @property
    def child_tags(self):
        return [x.tag.split('}')[1] for x in self]
    
    @property
    def text(self):
        return self._element.text

    def find(self, key:str):
        # 获取指定的子节点
        elm = self._element.find(self.tag_root + key)
        return self.__class__(elm)

    def findall(self, key:str):
        elms = self._element.findall(self.tag_root + key)
        return [self.__class__(x) for x in elms]
    
    def find_deep(self, tag:str):
        """
        查找指定标签名的元素，包括嵌套的子节点
        :param tag: 要查找的标签名
        :return: 包含所有匹配元素的列表
        """
        if not tag.startswith(self.tag_root):
            tag='{http://www.hmdb.ca}'+ tag
        result = []
        if self.tag == tag:
            result.append(self)
        for child in self:
            result.extend(child.find_deep(tag))
        return result
    
    def collect_phys_effect(self):
        # 找到当前记录中存方功能和疾病的root
        targets = self.find_deep('root')
        for it in targets:
            if it.get('term') == 'Physiological effect':
                result = it.find_deep('descendant')
                return [x.get('term') for x in result if x.get('type')=='child']  

    def collect_disease(self):
        data = []
        for node in tqdm(self, desc='processing', ncols=100):
            item = {}
            item['name'] = node.get('name')
            item['smiles'] = node.get('smiles')
            item['hmdb_id'] = node.get('accession')
            item['disease'] = [disease.get('name') for disease in node.find_deep('disease')]
            data.append(item)
        df = pd.DataFrame(data)
        df['n_ft'] = df['disease'].apply(len)
        print('converting ...')
        df = df[df['disease'].notnull()].explode('disease')
        return df[['hmdb_id', 'smiles', 'disease']].groupby('disease')['hmdb_id'].apply(list).reset_index()
    
    def collect_cmpd(self, fields = ['name', 'smiles', 'inchikey'], dropna= True):
        data = []        
        for node in tqdm(self, desc='processing', ncols = 100):
            it = {}
            it['hmdb_id'] = node.get('accession')
            for field in fields:
                it[field] = node.get(field)
            data.append(it)
        df = pd.DataFrame(data)
        if dropna:
            df = df[df['smiles'].notnull()]
        return df
    
    def get(self, tag):
        '''
        获得直接、特定子节点的值,只获得一个
        '''
        return self.find(tag).text
    
    def getall(self, tag):
        result = self.findall(tag)
        return [it.text for it in result]
    
    def base_df(self,
                fields:list = ['accession', 'name', 'smiles', 
                               'inchikey', 'cas_registry_number', 
                               'monisotopic_molecular_weight'],
                append:list = None) -> pd.DataFrame:
        if append is not None:
            fields += append
        result = []        
        for it in tqdm(self, desc='Processing', ncols=100):
            mol = {}
            for field in fields:
                mol[field] = it.get(field)
            result.append(mol)
        df = pd.DataFrame(result)
        df['monisotopic_molecular_weight'] = df['monisotopic_molecular_weight'].astype(float)
        return df     

class Hmdb():
    # cmpd = pd.read_pickle('../data/hmdb/hmdb_kegg_cmpd.pkl')

    @classmethod
    def load_disease(cls):
        disease = pd.read_pickle('../data/hmdb/disease.pkl')
        return Enrichment(disease)

    @classmethod
    def to_hmdb_id(cls, ids,
                on='kegg',
                dropna = False):
        '''
        Get hmdb_id with the list of on fields as index
        ids, list of id values.
        on, searching field.
        dropna, if to drop nan values in the result list
        '''
        if on =='kegg':
            on ='kid'
        if isinstance(ids, str):
            hmdb_id = cls.cmpd[cls.cmpd[on] == ids]['hmdb_id'].values
            if len(hmdb_id) == 0:
                return None
            elif len(hmdb_id) ==1:
                return hmdb_id[0]
            else:
                raise ValueError('More than one item was found!')
        else:
            if dropna:                
                return cls.cmpd[cls.cmpd[on].isin(ids)]['hmdb_id'].tolist()
            else:
                ids = pd.DataFrame({on: ids})
                result = pd.merge(ids, cls.cmpd, how='left', on='kid')
                return result['hmdb_id'].tolist()

