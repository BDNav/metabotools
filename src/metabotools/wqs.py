
import pandas as pd
from IPython import embed

def prep_for_wqs(infile, case, control, zscores=True):
    '''
    calculates zlogs and changes case/control to binary
    
    '''
    from correlation_analysis.correlation_analysis import get_zscores
    data = pd.read_csv(infile, index_col=0)
    
    if zscores==True:
        res = get_zscores(data, case, control)

        data2 = res['case'].append(res['control'])

        data2['Group'] = data['Group']

        data2 = data2[['Group'] + data2.columns[:-1].tolist()]

        
    else:
        data2 = data
    
    data2['Group'] = data2['Group'].replace(case,1)
    data2['Group'] = data2['Group'].replace(control,0)
    
    return data2

def run_wqs(data, num_bootstraps=10):
    '''
    takes in a transformed dataframe outputs weights from WQS
    
    '''
    import os
    transformed_file = 'tmp_for_wqrs.csv'
    
    colkey = dict(zip(data.columns[1:], range(0,len(data.columns[1:])))) # need to rename cols for input to R
    
    data.rename(columns=colkey, inplace=True)
    
    data.to_csv(transformed_file)

    R_script_dir = (os.sep).join(os.path.abspath(__file__).split(os.sep)[:-1])+os.sep+'R_scripts'
    
    cmd = 'Rscript.exe %s/WQS_regression.R "%s" %i'%(R_script_dir, transformed_file, num_bootstraps) # +1 because R index start at 1
    print(cmd)
    os.system(cmd)
    import time
    time.sleep(2)
    
    weights = pd.read_csv('weights.csv', index_col=0)
        
    weights['name'] = weights['mix_name'].str.replace('X','')
    weights['name'] = pd.to_numeric(weights['name'])
    
    colnums = pd.DataFrame({'colnum':colkey}).reset_index()
    
    out = pd.merge(weights, colnums, left_on='name', right_on='colnum')
    out.rename(columns={'index':'MRM Name','mean_weight':'WQS Weight'}, inplace=True)
    
    out['WQS Rank'] = range(1,len(out)+1)
    
    return out[['MRM Name','WQS Weight','WQS Rank']]
    
    
    
if __name__=='__main__':

    
    
    

    infile = 'D:\\Dropbox (Personal)\\naviauxlab_informatics\\Suicide_Study\\WQS\\wqs_no_na.csv'
    
    
    
    
    
    res = run_wqs(infile)

    embed()
    