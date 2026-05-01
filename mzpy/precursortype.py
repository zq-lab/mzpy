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

# from .mz import centroid




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

