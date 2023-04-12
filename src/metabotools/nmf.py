
import pandas as pd
from IPython import embed


def run_nmf(transformed_file):
    '''
    takes in a transformed dataframe outputs weights from WQS
    
    '''
    import os
    
    R_script_dir = (os.sep).join(os.path.abspath(__file__).split(os.sep)[:-1])+os.sep+'R_scripts'
    
    cmd = 'Rscript.exe %s/NMF.R "%s"'%(R_script_dir, transformed_file) # +1 because R index start at 1
    print(cmd)
    os.system(cmd)

    
    
if __name__=='__main__':

    
    
    

    infile = 'D:\\Dropbox (Personal)\\naviauxlab_informatics\\Suicide_Study\\April2022\\Plasma_M_suicide_Sho_rap_FGF+GDF_20190926_N=93_PPCA.csv'
    
    
    
    
    
    res = run_nmf(infile)

    