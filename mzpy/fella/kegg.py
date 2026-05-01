"""
Classes for handling KEGG record from rest.kegg.jp
' KEGG is a rule-based storage database. 
' KEGG IDs follow strict naming conventions. 
' This code package is built on the basis of KEGG rules.
' Glycon is not included !! 
' Glycon is not syitable for general analysis of metabolites.

# terms, conception or abbrevations:
    : kid, kegg identifier
    : kids, kegg identifiers
    : category, kegg class, kegg node type, including organism, pathway, module, reaction, enzyme and compound
    : dm, differential metabolites
    : node, kegg object,such as R00010, M00020
    : link, the connection between two kegg object
    : from, stat point of a link
    : to, end point of a link
    : weight, distance between two nodes
    : enr, enrichment

# code conventions
    : kegg.zip, kegg node txt downloaded from rest.kegg.jp
    : kegg_metabolites.csv, all metabolites list involved in kegg metabolic metwork.
    : kegg_networks.h5, contains basic dataframes for enrichment and igraph:
        : enr_base, a basic table for enrichment
        : links, a from-to table for igraph edges covering the whole network of kegg

# note
    : I once tried to download all the text of KEGG nodes to a local folder, 
        but it was easily rejected by KEGG for access. This is not a good solution.
    : For storing data, HDF5 is much better than SQLite
    : Module, reaction, enzyme, and compound are not species-specific.
    : The pathway represented by the map number (such as map00010) is not species-specific
        in terms of its function.
    : Pathways are classified by species, and the functions of pathways in different species are the same.
        But the modules and compounds they contain may differ.

# code reference
    : FELLA, an R package shared in Biocnductor community
"""
from abc import ABC, abstractmethod
import ast
from concurrent import futures
from datetime import datetime
import igraph as ig
from matplotlib import cm
import os
import pandas as pd
import random
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import re
import requests
import sqlite3
import time
from tqdm import tqdm
import urllib3

headers = {  
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',  
    'Accept-Language': 'en-US,en;q=0.9',  
    'Accept-Encoding': 'gzip, deflate, br',  
    'Connection': 'keep-alive',  
} 


def get(kid:str) -> str:
    # if a kid ends with '.mol', then download the mol text file.
    if kid.endswith('.mol'):
        url = 'https://www.kegg.jp/entry/-f+m+'
        url = url + kid.split('.')[0]
        response = requests.get(url)  
        if response.status_code == 200:
            return Chem.MolFromMolBlock(response.text) 
        else:
            raise ValueError(f'{response.status_code}:: response error, please check kid value')         
    else:
        url = f'https://rest.kegg.jp/get/{kid}'
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            return f'{response.status_code}:: response error'


def ls(category:str):
    '''
    :param 
        category, kegg category if you want list organism pathway,
                    such as "dre", the category should be set as 'pathway/dre'
    '''

    if category in Node.pattern.keys():
        df = pd.read_csv('https://rest.kegg.jp/list/'+category,
                         header=None,
                         names=['kid', 'name'],
                         sep='\t')
        return df
    elif  bool(re.match(r'^[a-z]{3,4}$', category)):
        '''为指定物种的缩写时'''
        df = pd.read_csv('https://rest.kegg.jp/list/pathway/'+category,
                         header=None,
                         names=['kid', 'name'],
                         sep='\t')
        return df
    elif category == 'organism':
        df = pd.read_csv('https://rest.kegg.jp/list/organism',
                         header=None,
                         names=['kid', 'name', 'latin_name', 'taxonomy'],
                         sep='\t')
        return df        
    else:
        raise ValueError('unrecognized category value!')     


def create(txt):
    '''
    auto-initialize sub class of kegg
    '''
    entry = txt.split('\n', 1)[0]
    entry = entry.split()
    if len(entry) > 2:  
        category = entry[-1]
        category = category.capitalize()   
        return eval(category)(txt)
    else:
        raise ValueError('Unrecoginzed kegg id, please check it.')



class LocalDb():
    def __init__(self, db_file:str):
        self.conn = sqlite3.connect(db_file)
        self.cursor = self.conn.cursor()
        self.db_file = os.path.abspath(db_file)
    
    def __del__(self):
        self.conn.close()

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__}: {self.db_file}>'
    
    def close(self):
        self.conn.close()

    ### reading
    def get(self, kid):
        '''
        get txt or molecule object according to the kid
        '''
        if kid.endswith('.mol'):
            sql = f'SELECT smiles FROM molecule WHERE kid = "{kid[:-4]}"'
            result = self.cursor.execute(sql).fetchall()
            if len(result) > 0:
                smiles = result[0][0]
                return Chem.MolFromSmiles(smiles)
            else:
                return None
        else:
            sql = f'SELECT txt FROM item WHERE kid = "{kid}"'
            result = self.cursor.execute(sql).fetchall()
            if len(result) > 0:
                return result[0][0]
            else:
                return None

    def get_cmpd_by_org(self, org:str):
        sql = f'SELECT compounds FROM features WHERE kid LIKE "{org}%"'
        df = pd.read_sql(sql, self.conn)
        df['compounds'] = df['compounds'].apply(ast.literal_eval)
        kids = df['compounds'].explode().tolist()
        return set(kids)

    def ls(self, category=None, org:str=None):
        '''
        category, can be a string or a list.
        '''      
        sql =  f'SELECT kid FROM item '
        if isinstance(category, str):
            sql += f'WHERE category = "{category}"'
        elif isinstance(category, list):
            sql += 'WHERE category IN '
            sql += str(category).replace('[', '(').replace(']', ')')
        result = self.cursor.execute(sql).fetchall()
        kids = [kid for sub in result for kid in sub]
        if (org is not None):     
            if category == 'pathway':      
                kids = [kid for kid in kids if kid.startswith(org)]
            else:
                raise ValueError('If you want get org pathway entries, category must set to "pathway"')
        return kids
    
    def load(self, tbl_name):
        df = pd.read(f'SELECT * FRoM {tbl_name}', self.conn)
        return df
    
    def load_ft(self, org, gene_included=False):
        '''
        应该逐层抓取数据
            先抓物种的通路，根据module信息，再抓module，然后是reaction，然后是enzyme
            所有物种的通路里，包含了物种全部的代谢物
            根据
        '''
        if gene_included:
            df = self.load('ft_cmpd_gene')
        else:
            df = self.load('ft_cmpd')
        df = df[df['category'].isin(['module', 'reaction', 'enzyme']) |
                df['kid'].str.startswith(org)]
        df['feature'] = df['feature'].apply(ast.literal_eval)
        return df

    ### save data
    def build_ft(self, org, gid=False, limit=0):
        '''
        构建用于富集分析的特征表
        '''
        ## 先处理 pathway
        print('processing pathway data...')
        sql = f"SELECT kid, txt FROM item WHERE kid LIKE '{org}%'"
        if limit > 0:
            sql = sql + f' LIMIT {limit}'
        self.cursor.execute(sql)
        result = []
        modules = []

        for row in self.cursor.fetchall():
            kid, txt = row
            print(f'{kid} ......', end='\r', flush=True)
            it = create(txt)
            result.append(it.exft(gid=gid))
            if hasattr(self, 'module'):
                modules = modules + list(it.module)

        ## 处理module
        print('\nprocessing module data ...')
        sql = "SELECT kid, txt FROM item WHERE category = 'module'"
        self.cursor.execute(sql)
        reactions = []
        for row in self.cursor.fetchall():
            kid, txt = row
            print(f'{kid} ......', end='\r', flush=True)
            if kid in modules:
                it = create(txt)
                result.append(it.exft(gid=gid))
                if hasattr(it, 'reaction'):
                    reactions = reactions + list(it.reaction)

        ## 处理reaction
        print('\nprocessing Reaction data ...')
        sql = "\nSELECT kid, txt FROM item WHERE category = 'reaction'"
        self.cursor.execute(sql)
        for row in self.cursor.fetchall():
            kid, txt = row
            print(f'{kid} ......', end='\r', flush=True)
            if kid in reactions:
                it = create(txt)
                result.append(it.exft(gid=gid))

        df = pd.DataFrame(result)
        df['n_feature'] = df['feature'].apply(len)
        return df[df['n_feature']>0]


    def check_txt_completeness(self):
        sql = f'SELECT kid FROM item WHERE txt NOT LIKE "%///%"'
        result = self.cursor.execute(sql).fetchall()
        return [kid for sub in result for kid in sub]

    def check_missing_items_in_category(self, category):
        '''Check if all data is downloaded and saved according to the kegg category.
        : return Return a list of difference kegg identifier'''
        full_list = ls(category)
        saved = self.ls(category)
        return (set(full_list) - set(saved))
            
    def download(self,
                 kids,
                 dtype:str='item', # ('item', 'mol')
                 n_thread = 1,
                 delayed:bool=False,
                 chunksize:int = 50) -> None:
        ''' 
        - The downloaded data is saved to the zip file specified by current object (self).
        - Make sure that write permissions are turned on for the object
        - if a kid ends with '.mol', then download the mol text file.
        param:
            kids, a list of kegg id for downloading
            dtype, item, general information txt; molecule, mol txt for molecules
                    it is also the strorage table name.
        '''

        cur_date = str(datetime.now().date())
        i = 0
        n_thread = min(20, len(kids), n_thread)
        with futures.ThreadPoolExecutor(max_workers=n_thread) as executor:
            i += 1
            random.shuffle(kids) # 该函数就地打乱kids的顺序，返回None
            for kid in tqdm(kids, desc='Processing', ncols=100):
                if delayed:
                    delay = random.uniform(0.5, 2) # randomly generate number from 0.5 to 2
                    time.sleep(delay)
                if dtype == 'mol':
                    kid = kid + '.mol'
                future = executor.submit(self.rest.get, kid)
                try:
                    txt = future.result()
                except Exception:
                    pass 
                else:
                    if dtype == 'mol':
                        cate = 'compound'
                    else:
                        cate = Node.identify(kid)
                        txt = txt.replace('"', "'")
                    sql = f'''INSERT OR REPLACE INTO {dtype} 
                        (kid, txt, update_on, category)
                        VALUES ("{kid[:6]}", "{txt}", "{cur_date}", "{cate}")'''
                    self.cursor.execute(sql)           
            if i % chunksize == 0:
                self.conn.commit()
            self.conn.commit()

    def download_by_category(self, category, n_thread=10):
        kids = self.rest.ls(category)
        self.download(kids=kids, n_thread=n_thread)

    def load(self, tbl_name):
        return pd.read_sql(f'SELECT * FROM {tbl_name}', self.conn)
    
    def append_df(self, df, sql_table_name, duplicated_on:list = None):
        '''append df to local database and de_duplicated'''
        try:      
            df2 = self.load(sql_table_name=sql_table_name)
        except:
            df2 = None
        data = pd.concat([df2, df]) # 即使df2为None，也可以合并
        if duplicated_on:
            data = data.drop_duplicates(subset=duplicated_on)
        n = data.to_sql(sql_table_name, self.conn, if_exists='replace', index=False)
        print(f'{df.shape[0]} items are incorprated into {sql_table_name}, {n} items in total.')
    
    def build_links(self, kids:list, sql_table_name:str = None):
        '''
        generate links table according to provided kids list.
            If sql_table_name is provided, it is stored in the local sqlite3 database.
        return: a dataframe containing 3 columns: from_kid, to_kid and category of to_kid.            
        '''
        data = []
        for kid in tqdm(kids, desc='Processing', ncols=100):
            xx = Node.create(kid)
            data += xx.to_links()
        df = pd.DataFrame(data, columns=['from_kid', 'to_kid', 'category'])
        df = df.drop_duplicates(subset=['from_kid', 'to_kid'])
        if (sql_table_name is not None) & (sql_table_name != ''):
            self.append_df(df, sql_table_name = sql_table_name, duplicated_on=['from_kid', 'to_kid'])
        return df
    
    def build_links_batch(self,
        category = ['pathway', 'module', 'enzyme', 'reaction', 'ko', 'disease'],
        org:str = None,
        sql_table_name:str = None):
        kids = self.ls(category = category, org = org)
        if (kids is not None) & (len(kids) > 0): 
            return self.build_links(kids=kids, sql_table_name=sql_table_name)
        else:
            print('no proper kegg id was found.')
    
    # def build_features(self, kids, sql_table_name:str = None):
    #     collection = []
    #     for kid in tqdm(kids, desc='Processing', ncols=100):
    #         xx = Node.create(kid)
    #         if 'compound' in xx:
    #             data = {}
    #             data['kid'] = kid
    #             data['category'] = xx.category
    #             if 'compound' in xx:
    #                 data['compounds'] = str(set(xx.compound))
    #             data['n_compounds'] = len(xx.compound)
    #             collection.append(data)
    #     df = pd.DataFrame(collection)
    #     if (sql_table_name is not None) & (sql_table_name != ''):
    #         self.append_df(df, sql_table_name=sql_table_name, duplicated_on=['kid'])
    #     return df
    
    # def build_features_batch(self,
    #     category = ['pathway', 'module', 'enzyme', 'reaction', 'ko'],
    #     org:str = None,
    #     sql_table_name:str = None):
    #     kids = self.ls(category = category, org = org)
    #     if (kids is not None) & (len(kids) > 0): 
    #         return self.build_features(kids=kids, sql_table_name=sql_table_name)
    #     else:
    #         print('no proper kegg id was found.')
    
    def load_ko_enrich_base(self, sql_table_name = 'features'):
        sql = 'SELECT * FROM features WHERE category = "ko"'
        df = pd.read_sql(sql, self.conn)
        df['compounds'] = df['compounds'].apply(ast.literal_eval)
        return df

    def load_pathway_enrich_base(self, org:str = None, sql_table_name = 'fefatures'):
        if org is None:
            sql = 'SELECT * FROM features WHERE category IN ("enzyme", "module", "reaction", "pathway")'
        else:
            sql =  'SELECT * FROM features WHERE category IN ("enzyme", "module", "reaction") '
            sql += f'OR kid LIKE "{org}%"'
        df = pd.read_sql(sql, self.conn)
        df['compounds'] = df['compounds'].apply(ast.literal_eval)
        return df
    
    def load_links(self, targets:list,
                    expand:str = 'half'):
        '''
        Obtain link data according to the kids in targets.
        expand, extension level of link (tree) data. ('full', 'half', 'rigid')
            full, all the "to_kid" of node in targets,
                    including all child and all descendant nodes.
            half, all to_kid in target and direct child nodes.
            rigid, from_kid and to_kid are limited in targets. 
        
        '''
        df = self.load('links')
        if expand == 'full':
            def extract(target):
                targets = [kid for kid in target if Node.identify(kid)!='compound']
                if len(target) == 0:
                    return None
                else:
                    data = df[df['to_kid'].isin(targets)]
                    return pd.concat([data, extract(data['from_kid'])])
            return extract(targets)
        elif expand == 'half':
            return df[df['to_kid'].isin(targets)]
        elif expand == 'rigid':
            return df[(df['from_kid'].isin(targets)) & (df['to_kid'].isin(targets))]
        else:
            raise ValueError('expand should be one of ("full", "half", "rigid").')

    def load_org_links(self, org, sql_table_name = 'links'):
        sql = 'SELECT * FROM links WHERE (category IN ("module", "reaction", "enzyme"))'
        sql += f'OR to_kid LIKE "{org}%"'
        df = pd.read_sql(sql, self.conn)
        return df
    


class Node(ABC):
    '''base class for KEGG item txt parase'''
    pattern = {
            # 'organism': r'[a-z]{3,4}',
               'pathway' : r'[a-z]{3,4}\d{5}',
               'module'  : r'M\d{5}',
               'enzyme'  : r'\d+\.[\d-]+\.[\d-]+\.[\d-]+',
               'reaction': r'R\d{5}',             
               'compound': r'C\d{5}',
               'ko'      : r'K\d{5}',
               'disease' : r'H\d{5}'}

    @property
    def categories(self):
        # property仅用于访问对象属性，不能用于访问类属性，
        return list(self.pattern.keys())

    @classmethod
    def identify(cls, kid):
        # identify the category based on kid
        for key in cls.pattern:
            if re.fullmatch(cls.pattern[key], kid):
                return key    
       
    def __init__(self, txt):
            self.txt = txt
            self._parse()
                
    def __contains__(self, item):
        return item in self.__dict__
    
    def __getitem__(self, key):
        return self.__dict__[key]
 
    def _parse(self):
        '''把kegg文本拆分为键值对数组'''
        tags = re.findall(r'(\A\w+ )', self.txt) + \
               re.findall(r'(\n\w+ )', self.txt)
        tags.append('///') # 全部记录结束符
        for i in range(len(tags)-1):
            pattern = re.compile(fr"{tags[i]}(.*?){tags[i+1]}",  re.DOTALL)
            result = pattern.search(self.txt)
            key = tags[i].strip().lower()  # 属性名（键）全部小写
            if result:
                self.__dict__[key] = result.group(1).strip()
            else:
                self.__dict__[key] = None

        if 'entry' in self.__dict__:
            entry = self.entry.split()
            self.kid = entry[0]
            self.category = entry[1]  
        if 'name' in self.__dict__:
            self.name = self.name.split(';', 1)[0]
        if 'gene' in self.__dict__:
            lines = self.gene.strip().split('\n')  # 按行分割  
            lines = [line.strip() for line in lines]
            gene_id = [re.match(r'^\d+', line).group(0) for line in lines if re.match(r'^\d+', line)]
            self.gene_id = [id for id in gene_id if (id is not None) and (id != '')]
            # 不是所有的基因都有symbol，所以以gene id作为依据
            # 例如dre04060中的基因：100329520、100332902、101887111都没有symbol

        # 如果包含其它的关联节点，则提取节点kid
        for key in self.pattern: 
            if key in self.__dict__:
                self.__dict__[key] = set(re.findall(self.pattern[key],
                                                self.__dict__[key]))

    def get(self, key, default_value=None):
        return self.__dict__.get(key, default_value) 
    
    def exft(self, gid=False):
        '''
        export freature, feature can be metabolite ID (kid) or gene id (gid)
        '''
        features = []
        if hasattr(self, 'compound'):
            features = features + list(self.compound)
        if gid and hasattr(self, 'gene_id'):
            features = features + list(self.gene_id)
        features = [id for id in features if (id is not None) and (id != '')]
        features = set(features) 
        return {'kid': self.kid, 'category':self.category, 'feature':features}    

    @abstractmethod
    def to_links(self, category:str) -> list:
        '''
        create the links from the specific category to self
        Return:
            a 2D list with 3 column: from, to and weight. Such as:
            [['C00078', 'M00010', 'module'],
             ['C00079', 'M00010', 'module']]        
        '''    
        links = []
        if category in self:
            for kid in self.__dict__[category]:
                links.append([kid, self.kid, self.category])
        return links
    
   
    def __repr__(self) -> str:
        return f'<KEGG {self.category}: {self.kid}, {self.name}>'    


# class Organism(Node):
#     # Essentially, it is a calss for processing pathways of a species-specific collection.
#     def __init__(self, org):
#         self.kid = None
#         self.name = org
#         self.kids = self.rest.ls(org)
    
#     def to_links(self) -> None:
#     # 屏蔽父类, 这类节点不应该有kegg网络
#         return []


class Pathway(Node):
    def __init__(self, txt):
        super().__init__(txt)
        if 'gene' in self:
            enzyme = re.findall(r'\[EC\:(.*?)\]', self.gene)
            enzyme = [item for row in enzyme for item in row.split(' ')]
            enzyme = set(enzyme) - set([''])
            self.enzyme = enzyme

    def to_links(self) -> list:
        links = []
        # 必须在以下四种分支中四选一，否则网络太复杂
        # 有的通路没有module，比如：cel00220
        if 'module' in self:
            links += super().to_links('module')
        elif 'enzyme' in self:
            links += super().to_links('enzyme')
        elif 'reaction' in self:
            links += super().to_links('reaction')
        elif 'compound' in self:
            links += super().to_links('compound')
        return links


class Module(Node):
    def __init__(self, txt):
        super().__init__(txt)   
        enzyme = re.findall(r'\[EC\:(.*?)\]', self.orthology)
        enzyme = [item for row in enzyme for item in row.split(' ')]
        enzyme = set(enzyme) - set([''])
        self.enzyme = list(enzyme)

    def to_links(self) -> list:
        links = super().to_links('enzyme')
        links += super().to_links('reaction')
        return links


class Enzyme(Node):
    def __init__(self, txt):
        super().__init__(txt) 
        if 'compound' not in self:
            if 'substrate' in self:
                self.compound = re.findall(self.pattern['compound'], self.substrate)
            if 'product' in self:
                self.compound += re.findall(self.pattern['compound'], self.product)
    
    def find_genes(self, org):
        '''
        Extracting gene information from the genes text.
        The gene ID in KEGG like hsa: 144811.
        the gene id is also NCBI-GeneID.
        '''
        if 'genes' not in self:
            return None

        org = org.upper()
        ss = re.findall(rf'^\s*{org}: (.*)$', self.genes, re.MULTILINE)
        if (ss is None) or (len(ss) ==0):
            return [{'enzyme' : self.kid,
                     'gene_id': '',
                     'symbl'  :''
            }]
        else:
            data = []
            for s in ss[0].split(' '):
                g = {}
                g['enzyme'] = self.kid
                g['gene_id'] = s.split('(')[0].strip()
                symbls = re.findall(r'\((.*?)\)', s)
                g['symbl'] = symbls[0] if len(symbls) > 0 else ""
                data.append(g)
            return data

    def to_links(self) -> list:
        return super().to_links('reaction')

    @classmethod
    def get_assoc_genes_seq(cls, org:str, kids:list):
        '''
        According to the given kids of enzymes,
            obtain corresponding genes information and sequence.
        '''
        data = []
        for id in kids:
            enz = cls(id)
            og = enz.find_genes(org)
            if og is not None:
                data = data + enz.find_genes(org)
        for it in data: #给data的每个条目增加aaseq记录
            if ('gene_id' in it) and (it['gene_id'] != ''):
                gene = Gene(org, it['gene_id'])
                it['aaseq'] = gene.aaseq
            else:
                it['aaseq'] = ''
        return pd.DataFrame(data)


class Reaction(Node):
    def __init__(self, txt) -> None:
        super().__init__(txt)   
        self.compound = re.findall(self.pattern['compound'], self.equation)
    
    def to_links(self):
        return super().to_links('compound')
  

class Compound(Node):   
    def load_mole(self):
        molTxt = self.zipf.get_mol_block(self.kid)
        try:
            mole = Chem.MolFromMolBlock(molTxt)
            self.mole = mole
            return mole
        except:
            return None
    
    def to_links(self) -> list:
        return [] # 应该与其他兄弟类保持形同类型的输出

class Ko(Node): # 字母o是小写
    '''
    KEGG has a total of 26,410 records (2024-12-20). 
    Ther are 20,521 KO records not involved reactions, accounting for 77% in total.
    '''
    def __init__(self, txt):
        super().__init__(txt)  
        if 'reaction' in self:
            compound = []
            for kid in self.reaction:
                xx = super().create(kid)
                compound += list(xx.compound)
            self.compound = list(set(compound))

    def to_links(self):
         return super().to_links('reaction')
    
    def get_features(self):
        data = []
        if 'reaction' in self:
            for kid in self.reaction:
                xx = super().create(kid)
                data += xx.compound
        return data



class Disease(Node): 
    def __init__(self, txt):
        super().__init__(txt)
        if 'gene' in self:
            self.ko = re.findall(self.pattern['ko'], self.gene)
    def to_links(self):
        return super().to_links('ko')


class Gene(Node):
#     # the url (https://rest.kegg.jp/list/genes) does not work.
#     # The sequence information of gene can only be obtained online. 
#     # Since genes list cannot be obtained, so it cannot be downloaded locally.
#     # via url: https://rest.kegg.jp/get/hsa:144811
    def __init__(self, org, kid):
        # 只能在线获取数据
        kid = f'{org}:{kid}'
        super().__init__(kid, from_local=False)
    
    def to_links(self):
        return



class Net():
    '''
    Drawing a network diagram, which only needs two steps: 
        - Initialize objects from a from... to table 
        - Call the plot function 
    '''
    style  = {'vertex_label_size': 9,
              'vertex_label_dist': 0.5,
              'edge_color'       : 'black',
              'edge_curved'      : 0,
              'edge_width'       : 1.0}
    
    shapes = ('circle','square', 'triangle', 'diamond', 'ellipse', 'polygon',
              'star',  'image',  'text',     'arrow',   'no_shape')
    
    sizes = {'pathway' : 26,
             'module'  : 20,
             'enzyme'  : 12, 
             'reaction': 12, 
             'compound': 8,
             'masked'  : 3} # masked vs
    
    # palette = dict(zip(Node.pattern.keys(), cm.__dict__['Set2'].colors))
    
    def __init__(self, 
                 targets:list, 
                 expand = 'half',
                 directed:bool = False,
                 masked = True,
                 cate_size = None,
                 palette = None, 
                 *args, **kwargs):
        # 下一代增加mask参数，mask = True，则不在targets中的节点用灰色
        # 在富集的网络图中，不在富集结果的节点显示灰色，便于区分
        '''
        Filter the links list based on the supplied node list (to_kids)
            and initialize the igraph object.
        param:
            to_kids, a list of to_kid in link table
            expand_tree, 是否包含子节点的子节点
            shape: a list for kid or node id, which need marked as square.
        '''
        db = LocalDb()
        links = db.load_links(targets=targets, expand=expand)
        links = links.iloc[:,0:2].values.tolist()
        db.close()       
        self.ig = ig.Graph.TupleList(links, directed = directed, *args, **kwargs)
        # TupleList(links.itertuples(index=False), directed=True, weights=False)
        # 使用links.itertuples(index=False)似乎能够绘制更完整的网络图
        self.ig.vs['category'] = [Node.identify(it) for it in self.ig.vs['name']]  
        self.ig.vs['size'] = [self.sizes[it] for it in self.ig.vs['category']]      
        if masked:
            self.mask_vs(targets=targets, palette=palette)
    
    @property
    def vs(self):
        return self.ig.vs # 这样也可以通过self['---']增加属性
    
    def mask_vs(self, targets:list = None, palette:str = 'Set2'):
        '''
        Target nodes are set to color according to the provided palette by category, 
            and the remaining nodes are set to gray.
        param:
            targets, nodes for coloring
            masked, to mask colors to grey for nodes not in targets.
            cate_size: node size  which set along their catgory.
            palette: palette for nodes
        '''
        if 'category' not in self.ig.vs.attributes():
            raise ValueError('vs category of iGraph is not set. category must be set firstly')
        self.palette = dict(zip(Node.pattern.keys(), cm.__dict__[palette].colors))
        color = []
        labels = []
        if targets:
            for name, cate in zip(self.ig.vs['name'], self.ig.vs['category']):
                if name in targets:
                    color.append(self.palette[cate]) 
                    labels.append(name)
                else:
                    # color.append((0.35, 0.35, 0.35, 0.75)) # grey
                    # the 4th value of RGB is alpha
                    color.append(cm.__dict__[palette].colors[-1])
                    labels.append('')
            self.ig.vs['color'] = color
            self.ig.vs['label'] = labels
        else:
            self.ig.vs['color'] = [self.palette[cate] for cate in self.ig.vs['category']]
            self.ig.vs['label'] = self.ig.vs['name']
            

    def set_label(self, targets, fontsize):
        return 
           
    def __str__(self):
        return self.g.__str__()
    
    # def set_palette(self, theme='Set1'):
    #     if isinstance(theme, str):
    #         self.palette = cm.__dict__[theme].colors
    #     elif isinstance(theme, tuple) or isinstance(theme, list):
    #         self.palette = theme #RGB list would be OK
    #     else:
    #         raise TypeError('Only str, flatten tuple or list can be acceptable!')
    
    def centrality(self):
        '''
        The most important enzyme is often around the most important module.
        '''
        df = pd.DataFrame({'name'      : self.ig.vs['name'],
                           'category'  : self.ig.vs['category'],
                           'centrality': self.ig.evcent()
            })
        df = df.sort_values(by='centrality', ascending=False)
        return df

    def plot(self, layout='fr', save_to=None):
        # 适合代谢通路的layout有fr和tree
        if 'label' not in self.vs.attribute_names():# display label
            self.vs['label'] = self.vs['name']
        return ig.plot(self.ig,
                        layout = self.ig.layout(layout),
                        target = save_to, # 输出的文件名
                        **self.style)
 



