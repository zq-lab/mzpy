##算法太繁琐，舍弃？
## 还是留着用于探索新型的precursor type？

'''
Function: To find appropriate MS1 peak with MSMS
Peaklist Files: exported from MS-Dial with version 5.5 or higher
How to find an appropriate MS1 peak?
    In general,  the cleavage only occured around heteroatom (alpha or beta position). 
    The smaller one will be lost as neutral fragment, and the bigger one will be charged as:
        [M-A]+-, [M-A-H]+- or [M-A+H]+-.
        [M-A+H]+-, which is equivalent to rearrangement, such as Maxwell rearrangement
    Two bonds can be broken in cyclic compounds to generated two fragments mostly.
        charged ions also will be [M-A]+-, [M-A-H]+-, or [M-A+H]+-.
    Be caution:
        The above does not cover all cases.
        The reliability needs to be tested according to a large amount of data.
        It is often difficult to find suitable types for molecules containing phosphate, 
            long chain molecules or steroids.

    Fragments can be divided into Three classes:
        1) the fragments contained in the molecule itself which we called "loss":
            such as H, OH, NH2, glycosyl, etc., which can be lost as neutral fragment
        2) adduct fragments from solvents, we called "adduct":
            such as H, H2O, Cl-, AcO-, etc
        3) loss or gain of an additional electron HCD: e

    Therefore, the determination of precursor type for organic molecules can be divided into 3 steps:
        loss a fragments form itself
        add an adduct from solvent
        add or loss an electron.
    
    Precursor Type:
        M + adduct_ion, M+H, M+Na
        M + adduct_ion + H(Na) or M + adduct_ion-H(Na)
        M - loss
        M - loss + H(Na)
        2M+H(Na)

'''

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors 
import re

from .ms import centroid

def is_pos(ionmode):
    return ionmode.lower() in ('+', 'p', 'pos', 'positive')

def is_neg(ionmode):
    return ionmode.lower() in ('-', 'n', 'neg', 'negative')

def chg(ionmode):
    '''
    judge ionmode
    return:
        'pos' or 'neg'
    '''
    if is_pos(ionmode):
        return 'pos'
    elif is_neg(ionmode):
        return 'neg'
    else:
        ValueError(f'unrecognized ionmode: {ionmode}')



data_folder = os.path.dirname(os.path.abspath(__file__))
mass = pd.read_csv(data_folder + '/mass.csv', comment='#')

def load_precursors(exact_mw, mode = None):
    '''
    获得不同类型母离子的类型和mz列表
    mode, None: all pos and neg ion type were calcd.
          pos, only calcd for pos ions
          neg, only calcd for neg ions
    '''
    if str(mode).lower() in ('+', 'p', 'pos', 'positive'):
        df = mass[mass['ionmode'] == 'Positive'].copy()
    elif str(mode).lower() in ('-', 'n', 'neg', 'negative'):
        df = mass[mass['ionmode'] == 'Negative'].copy()
    else:
        raise ValueError('unknown ionmode!')

    df['mz'] = df.apply(lambda row: \
                (row['n']*exact_mw + row['offset'])/row['charge'],
                axis=1)
    df['mz'] = df.apply(lambda row: (row['mz'] + 0.00055) \
                if row['ionmode'] == 'neg' else (row['mz'] - 0.00055),
                axis=1)
    df['mz'] = df['mz'].round(5)
    return df

def mz_matched(mz1, mz2, abs_thd = 0.003, rel_thd = 5):
    '''
    Determine whether the two mz values are within 
        the allowable margin of error. It is the main criteria 
        for MS1 Matching
    Smaller molecules are more suitable for absolute error measurements,
        while larger molecules are better assessed using relative error.
        Meeting either condition is sufficient.
    param:
        mz1, mz2: mz values
        abs_thd: Threshold (criteria) in absolute error, 0.003 Da
        rel_thd: Threshold (criteria) in relatvie error, 5 ppm
    return:
        True or False
    '''
    abs_err = abs(mz1 - mz2)
    rel_err = 1E6*abs_err/mz2
    return (abs_err < abs_thd) or (rel_err < rel_thd)


# 查找可断裂的化学键
def find_cleavage_bonds(mol, with_begin_end_idx = False):
    '''
    including alpha and beta bonds
    not inculding bonds of C-H, O-H, N-H etc.
    杂原子或不饱和碳原子alpha位单键
    '''
    bonds = [b for b in mol.GetBonds() if b.GetBondType() == \
                Chem.rdchem.BondType.SINGLE]
    targets = []
    for bond in bonds:
        begin = bond.GetBeginAtom()
        end   = bond.GetEndAtom()
        begin_idx = begin.GetIdx()
        end_idx = end.GetIdx()
        neighbor_atoms = list(begin.GetNeighbors()) + list(end.GetNeighbors())
        if begin.GetSymbol() != end.GetSymbol():
            targets.append((bond.GetIdx(), begin_idx, end_idx))
            continue
        elif begin.GetHybridization() != end.GetHybridization():
            if (begin.GetHybridization() == Chem.rdchem.HybridizationType.SP3) or\
                (end.GetHybridization() == Chem.rdchem.HybridizationType.SP3):
                targets.append((bond.GetIdx(), begin_idx, end_idx))
                continue
        for atom in neighbor_atoms:
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                targets.append((bond.GetIdx(), begin.GetIdx(), end.GetIdx()))
                break
            elif (atom.GetSymbol() != begin.GetSymbol()) or \
                (atom.GetSymbol() != end.GetSymbol()):
                targets.append((bond.GetIdx(), begin.GetIdx(), end.GetIdx()))
                break
    if with_begin_end_idx:
        return targets
    else:            
        return [idx for idx, _, _ in targets]

    


###===============================================
class PrecursorType():
    '''
    ## MoNA数据库中的precursor type或可提供更多的提示
    working like a calculator, which can only be used after initialization as an object.

    '''

    # e = 0.00055
    # H = 1.00783
    # H_mol = Chem.MolFromSmiles('[H]')

    # # mass does not contain the electrons
    # pos = {
    #     'H'    :  1.00783,
    #     'Na'   :  22.98977, 
    #     'K'    :  38.96371, 
    #     'NH4'  :  18.03437
    # }
    # neg = {
    #     'Cl'   :  34.96885, 
    #     'Br'   :  78.91834,
    #     'HCOO' :  44.99765,
    #     'AcO'  :  59.01330
    #     }
    # solvent = {
    #     'H2O'  :  18.01056, 
    #     'FA'   :  46.00548,
    #     'HAc'  :  60.02113,          
    #     'ACN'  :  41.02655,
    #     'TFA'  : 113.99286,
    #     'iPrOH':  60.05751,
    #     'DMSO' :  78.01394, 
    #     'MeOH' :  32.02621
    # }


    def display_idx(self, mol):
        '''
        show the mol structure with atom index
        '''
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return mol
    
    def get_charge(self, mol):
        '''
        get charge of a mol or fragment
        '''
        charges = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
        return sum(charges)

    def depart_complex(self, mol, ionmode):
        '''
        Separate organic salts and return ion information according to ionmode
        '''
        
        pass


 
    def generate_loss(self, mol, ionmode): 
        '''
        handling loss fragment
        '''
        ms1_types = {}
        recorded_loss_mf = []
        # 如果loss片段的mf相同，产生离子和衍生离子的mass也相同
        # 因此已经收集过的loss_mf，没必要再收集
        mass = Descriptors.ExactMolWt(mol)

        if chg(ionmode) == 'neg':
            ms1_types['[M-H]-'] = mass - self.H + self.e
            for solvent in self.solvent:
                ms1_types[f'[M+{solvent}-H]-'] = mass + self.solvent[solvent] -\
                        self.H + self.e

        sites = self.find_cleavage_bonds(mol)     
        for id in sites:
            fg = Chem.FragmentOnBonds(mol, [id], addDummies=True) 
                # 如果不加dummy原子会在断键处自动补氢原子
            try:
                fg = Chem.GetMolFrags(fg, asMols=True)
            except Exception as e:
                print(e)
                continue
            if len(fg) != 2:
                continue
            products = [{'mol': mol, 'mass':Descriptors.ExactMolWt(mol)} \
                for mol in fg]
            # 大片段作为产物，小片段作为丢失
            loss = min(products, key=lambda x: x['mass'])
            mf = rdMolDescriptors.CalcMolFormula(loss['mol'])
            if mf in recorded_loss_mf:
                continue
            else:
                recorded_loss_mf.append(mf)
            loss['mf'] = re.sub('\*\d*', '', mf)
            # 脱碎片后继续脱氢
            loss_H = {}
            loss_H['mol'] = Chem.CombineMols(loss['mol'], self.H_mol)
            mf_H = rdMolDescriptors.CalcMolFormula(loss_H['mol'])
            loss_H['mf'] = re.sub('\*\d*', '', mf_H)        
            if chg(ionmode) == 'pos':
                ms1_types[f'[M-{loss["mf"]}]+'] = mass - loss['mass'] - self.e
                ms1_types[f'[M-{loss_H["mf"]}]+'] = mass-loss['mass']-self.H-self.e
                for frag in self.pos:
                    ms1_types[f'[M-{loss["mf"]}+{frag}]+'] = \
                        mass - loss['mass'] + self.pos[frag] - self.e             
            elif chg(ionmode) == 'neg':
                ms1_types[f'[M-{loss["mf"]}]-'] = mass - loss['mass'] + self.e
                ms1_types[f'[M-{loss_H["mf"]}]-'] = \
                    mass - loss['mass'] - self.H + self.e
                for frag in self.neg:
                    ms1_types[f'[M-{loss["mf"]}+{frag}]+'] = \
                        mass - loss['mass'] + self.neg[frag] + self.e
        return ms1_types
    
    def generate_types(self, mol, ionmode):
        ms1_types = {}   
        if (self.get_charge(mol) > 0) and (chg(ionmode) == 'pos'):
            ms1_types['[M]+'] = Descriptors.ExactMolWt(mol)
            # 由于mol本身带着电荷，ExactMolWt在计算时已经考虑了电荷的因素
        elif (self.get_charge(mol) < 0) and (chg(ionmode) == 'neg'):
            ms1_types['[M]-'] = Descriptors.ExactMolWt(mol)
        elif self.get_charge(mol) == 0:       
            ms1_types.update(self.generate_adducts(mol, ionmode))
            ms1_types.update(self.generate_loss(mol, ionmode))            
        return ms1_types
        
    def is_matched(self, 
                    mass:float, 
                    mz:float,
                    delta:float = 7.0,
                    relative:bool = True):
        '''
        if error type is relative, delta was used in ppm
        '''
        if relative:
            error = 1.0E6 * abs(mass - mz) / mass
            return error < delta
        else:
            return (abs(mass - mz) * 1E3 ) < delta
        

    def identify(self, ms1_types:dict, mz:float):
        '''Determining the ion that matches the current mz according to 
                the given ions list.
        '''
        for t in ms1_types:
            if self.is_matched(mass = ms1_types[t], mz = mz):
                return t
        # if nothong found, return None

####### 为一个mzFram更新 precursortype
    def identify_types_df(self, df, mol, ionmode):
        mass = Descriptors.ExactMolWt(mol)
        if not all(col in df.columns for col in ['precursormz', 'intensity']):
            raise ValueError('Data frame must contains precursormz and intensity')
        # 获取TIC的峰尖的mz值
        ms1_centroid = centroid(df[['precursormz', 'intensity']].values)
        centroid_mz = [mz for mz, _ in ms1_centroid]

        df = df[df['Num Peaks'] > 0] 

        ms1_types =  self.generate_types(mol, ionmode=ionmode)
        df['precursortype'] = df['precursormz'].apply(lambda x: self.identify(ms1_types, x))

        for idx in df.index:
            precursortype = ''
            if pd.isnull(df.loc[idx, 'precursortype']):
                mz = df.loc[idx, 'precursormz']
                if  mz in centroid_mz:
                    if abs(mz - mass) < 5E-3:
                        precursortype = '[M]'
                    elif mz < mass:
                        precursortype = '[M-unk]'
                    else:
                        precursortype = '[M+unk]'

                    if chg(ionmode) == 'pos':
                        precursortype += '+'
                    else:
                        precursortype += '-'
                    df.loc[idx, 'precursortype'] = precursortype

        return df[df['precursortype'].notnull()]
