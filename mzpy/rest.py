import json
import pandas as pd
import re
import requests
from tqdm import tqdm
import urllib

# 设置用户代理（User-Agent）和其他头部信息
# 模拟浏览器访问
headers = {  
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',  
    'Accept-Language': 'en-US,en;q=0.9',  
    'Accept-Encoding': 'gzip, deflate, br',  
    'Connection': 'keep-alive',  
} 


def np_classfy(smiles):
    if (smiles is None) or (smiles == '') or (not isinstance(smiles, str)):
        return {'family': None,
        'superclass': None,
        'pathway': None}
    
    smiles = urllib.parse.quote(smiles)  
    url = f'https://npclassifier.gnps2.org/classify?smiles={smiles}'
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        info = json.loads(response.text)     
        family     = info['class_results'][0]      if len(info['class_results']) > 0      else ''
        superclass = info['superclass_results'][0] if len(info['superclass_results']) > 0 else ''
        pathway    = info['pathway_results'][0]    if len(info['pathway_results']) > 0    else ''
    elif response.status_code == 500:
        family = superclass = pathway  = 'unknown'
    else:
        family = superclass = pathway  = f'connection error: {response.status_code}'
    return {'family': family,
            'superclass': superclass,
            'pathway': pathway}


def np_classfy_df(df, smiles_on):
    np_class = []
    for idx in tqdm(df.index):
        smiles = df.loc[idx, smiles_on]
        class_info = np_classfy(smiles)
        class_info['smiles'] = smiles
        np_class.append(class_info)  
    
    return pd.DataFrame(np_class)  


class Compound:
    __slots__ = ['id', 'cas', 'name',
                 'inchikey', 'smiles', 'mf', 'mass',
                 'ontology']
    cas_pattern = r'\d{2,9}-\d{1,3}-\d'
    kid_pattern = r'C\d{5}'

    def __init__(self, id:str = None, cas:str = None, name:str = None,
                        inchikey:str = None, smiles:str = None, mf:str = None, mass:float = 0,
                        ontology:str = None):
        self.id         = id
        self.cas        = cas
        self.name       = name
        self.inchikey   = inchikey
        self.smiles     = smiles
        self.mf         = mf
        self.mass       = mass
        self.fill_ontology()


    @classmethod
    def create_from_pubchem(cls, field, key, show_url = False, structure_retrived=True):
        '''
        param:
            field, can be name, cid, cas or inchikey
            key, search term
            output
        统一数据的出入方式
        compound property请参考
        https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Compound-Property-Tables
        有些物质手工可以检索到但是rest检索不到, 即便是使用cid, rest依然检索不到, 原因不清楚, 比如下面的物质:
            QKCQVNOJOHBRGO-IRXDYDNUSA-N
            OIVAJJZZAYFRRZ-ONTIZHBOSA-N
            OTGHTIDLDTWPFJ-UHFFFAOYSA-N
            MWSJSUZPMVLSKC-UHFFFAOYSA-N
            RLHRLBNJCYNCSY-UHFFFAOYSA-N
            PZQCYAJOGCQXOA-UHFFFAOYSA-N
            UNBOSJFEZZJZLR-CCEZHUSRSA-N
            VFOJTRKXCVUSLX-FIXSFTCYSA-N
            UIWVAHSBXKQIOX-INIZCTEOSA-N
            SQVRNKJHWKZAKO-CPMYNFSBSA-N
            YABGGRNPJXBMFQ-ONTIZHBOSA-N
            有待查找原因优化代码
        '''
        if field not in ('cid', 'name', 'cas', 'inchikey'):
            raise ValueError("field must be one of ('cid', 'name', 'cas', 'inchikey')")  
        
        url_root = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/{}/json'
        key = str(key).replace(' ', '%20').replace('/', '%2F') 
        compound = cls()

        # 获取名称信息  
        url = url_root.format(field, key, 'synonyms')        
        if show_url:
            print(url)         
        response = requests.get(url, headers=headers)
        
        if response.status_code == 200:
            data = json.loads(response.text)      
            if 'InformationList' in data:
                data = data['InformationList']['Information']
                for it in data:
                    if 'Synonym' in it: 
                        # 有的inchikey检索后, pubchem返回多个cid值，但往往只有一个携带了Synonym
                        compound.name = it['Synonym'][0]
                        compound.id = str(it['CID'])
                        cas = [s for s in it['Synonym'] if re.match(cls.cas_pattern, s)]
                        compound.cas = ','.join(cas)
                        break


        #获取结构信息
        if compound.id and structure_retrived:
            url = url_root.format('cid', compound.id,
                                'property/Title,InChIKey,IsomericSMILES,MolecularFormula,MonoisotopicMass')
            if show_url:
                print(url)  
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                data = json.loads(response.text)
                if 'PropertyTable' in data:
                    values = data['PropertyTable']['Properties'][0]
                    compound.smiles   = values['IsomericSMILES']
                    compound.inchikey = values['InChIKey']
                    compound.mf       = values['MolecularFormula']
                    compound.mass     = values['MonoisotopicMass']
                    if 'Title' in values: #pubchem的可能有少数条目
                        compound.name = values['Title']
                response.close() # 可以先关闭连接，避免批量调用后爆内存
                compound.fill_ontology()

        return compound
    
    @classmethod
    def create_from_KEGG(cls, field, key, show_url=False):
        pass
    

 


            

        

