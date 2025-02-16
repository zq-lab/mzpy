import os
import pandas as pd


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