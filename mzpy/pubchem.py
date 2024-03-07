'''
Find compound information based on cid, inchikey, or name 
    - Can return a single value or a list as required 
    - Predefined an object (pubchem) that can be called directly 
'''



import urllib3

class PubChem():
    http = urllib3.PoolManager()
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/%s/%s/%s/TXT'
    # url = f'{url_root}/{by}/{key}/{field}/TXT'
    # abbrevations and standard terms used in the url
    field = {'cid'     : 'cids',
             'cids'    : 'cids',
             'name'    : 'Synonyms', 
             'smiles'  : 'property/IsomericSMILES',
             'inchi'   : 'property/InChI',
             'inchikey': 'property/InChIKey',
             'mf'      : 'property/MolecularFormula',
             'mw'      : 'property/MolecularWeight',
             'ms'      : 'property/MonoisotopicMass'}

    def search(self, by:str, key:str, to:str, first_only:bool=True, print_url:bool=False):
        '''
        - Searched according to cid, name or inchikey.
        - Others are not so good for searching'''
        if by in ('cid', 'name', 'inchikey'):
            full_url = self.url%(by, key, self.field[to])
            if print_url:
                print(full_url)
            return self.http_get(full_url, first_only=first_only)
        else:
            raise ValueError('search by must be "cid", "name" or "inchikey"!')

    def http_get(self, url:str, first_only:bool=True):
        if len(self.http.pools) > 10:
            self.http.clear()
        response = self.http.request("GET", url, timeout = 5)
        if response.status == 200:
            txt = response.data.decode().strip()
            data = txt.split('\n')
            if first_only and (len(data) > 0):
                return data[0]
            else:
                return data
        elif response.status == 404:
            return 'No record in pubChem'
        else:
            return 'Connection Error: %d'%response.status
