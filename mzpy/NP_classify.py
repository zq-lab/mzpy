# 获取天然产物分类信息
#！NPClassfier官方服务器不支持 futures多线程获取数据

from concurrent import futures
import urllib3
import json

http = urllib3.PoolManager()
url = 'https://npclassifier.ucsd.edu/classify?smiles='

def npClassify(smiles):
    smiles = smiles.replace('%', '%25').replace('+', '%2B').replace('#', '%23').replace('/','%2F')
    try:
        response = http.request('GET', url + smiles, timeout = 5)
    except:
        return {'status': -1,'class': '','superclass': '', 'pathway': ''}
    if response.status == 200:
        inf = json.loads(response.data.decode('utf-8'))     
        npClass      = inf['class_results'][0]      if len(inf['class_results']) > 0      else ''
        npSuperclass = inf['superclass_results'][0] if len(inf['superclass_results']) > 0 else ''
        npPathway    = inf['pathway_results'][0]    if len(inf['pathway_results']) > 0    else ''
        return {'status': 200, 'class':npClass,'superclass':npSuperclass, 'pathway':npPathway}
    else:
        return {'status': response.status, 'class': '','superclass': '', 'pathway': ''}  
    
def npClassifyDf(df, smiles_column, renew=False, n_thread=10):
    '''
    :param data, the working DataFrame
    :type data, DataFrame
    :param smiles_column, the colomn name of smiles string
    :type smile_column, str
    :param renew: if True, renew the whole table, 
                    if False, only working on those rows with status is not 200
    :type, bool, True or False
    '''
    # if df does not contain the column status, it will renew also 
    if renew or ('status' not in df.columns):
        df['status'] = 0

    i = 0
    with futures.ThreadPoolExecutor(max_workers=n_thread) as executor:
        for idx in df[df['status']!=200].index:
            i += 1
            kid = df.loc[idx, 'kid']
            print('%s is being processed ... ... %.1f %%'%(kid, 100*i/df.shape[0]),
                end='\r', flush=True)
            # info = npClassfy(df.loc[idx, smiles_column])
            future = executor.submit(npClassify, df.loc[idx, smiles_column])
            info = future.result()
            df.loc[idx,'status']      = info.get('status', '')
            df.loc[idx, 'pathway']    = info.get('pathway', '')
            df.loc[idx, 'superclass'] = info.get('superclass', '')
            df.loc[idx, 'class']      = info.get('class', '')
    return df