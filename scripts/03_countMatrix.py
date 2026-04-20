import os
import sys
import pandas as pd
from collections import defaultdict

def count_importer(count_file):
    df = pd.read_csv(count_file,
                       sep='\t',
                       comment='#')
    return df.set_index('Geneid')

def count_merger(count_files):
    '''
    Takes a list of count files generated for the same SAF
    files but different BAMs

    Parameter:
    ----------
    count_files : list
    A list of paths to .count files

    Returns:
    --------
    A Pandas DataFrame merging counts from each file
    based on Geneid. Because Geneid is used as the index,
    it is essential that all files were generated from the 
    same SAF.
    '''
    df_ls = [count_importer(f) for f in count_files]
    df_all = df_ls[0].copy()
    df_all = df_all.rename(columns={df_all.columns[-1]:\
                                    os.path.basename(df_all.columns[-1]).split('_')[0]})
    for df in df_ls:
        count_col = df.columns[-1]
        count_col = os.path.basename(count_col).split('_')[0]
        df_all[count_col] = df.iloc[:,-1]

    df_all['peakRegions'] = df_all['Chr'] + ':' + df_all['Start'].astype(str) + \
                            '-' + df_all['End'].astype(str)
    
    df_all = df_all.set_index('peakRegions',drop=True)
    df_all = df_all.drop(['Chr','Start','End','Strand','Length'],axis=1)
    df_all = df_all[sorted(df_all.columns)]
    
    return df_all

if __name__ == '__main__':
    args = sys.argv
    if len(args) >= 2:
        count_files = args[1:-1]
        output = args[-1]
        count_merger(count_files).to_csv(output)
    else:
        print("Usage: 03_countMatrix.py <count_files> <output>", file=sys.stderr)
        sys.exit(1)