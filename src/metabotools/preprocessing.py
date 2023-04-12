

# scripts written 6/17/2021 to impute missing data before input to downstream analytics (e.g. Table1, etc)


from IPython import embed


import pandas as pd

def main():

    # data = pd.read_csv('1_Plasma_exercise_PreR2_M+F_exploratory_4v8.csv') # test file
    
    input_file = 'Davis MIND Autism Pre-Vax_M+F_Met+Exp.csv'
    data = pd.read_csv(input_file)
    
    data = drop_too_many_na_cols(data, percent=100) # percent = drop features with > xx missing values  e.g. 100 -> keep all cols that are not ALL null
    
    # fill by mean of the column:
    '''
    data2 = fill_na_col_mean(data)
    '''
    
    # fill using ppca
    # '''
    # data2 = fill_na_ppca(data)
    # '''
    # embed()
    
    # fill using MetaboAnalystR:
    # requires to activate metaboanalyst2 via conda
    '''
    import os
    from shutil import move
    
    case = data[data['Group']=='ASD']
    # case.loc[case.index[:5],'Group'] = 'test'
    control = data[data['Group']=='TD']
    # control.loc[control.index[:5],'Group'] = 'test'
    
    case.to_csv('case.csv', index=False)
    control.to_csv('control.csv', index=False)
    
    cmd = 'Rscript.exe C:/Users/Jon/Documents/BITBUCKET/metabotools/R_scripts/ppca_impute.R "%s"'%'case.csv'
    print(cmd)
    os.system(cmd)
    move('data_processed.csv','case_processed.csv')
    
    cmd = 'Rscript.exe C:/Users/Jon/Documents/BITBUCKET/metabotools/R_scripts/ppca_impute.R "%s"'%'control.csv'
    print(cmd)
    os.system(cmd)
    move('data_processed.csv','control_processed.csv')
    
    case_transformed = pd.read_csv('case_processed.csv', index_col=0).transpose()
    case_transformed['Label'] = 'ASD'
    control_transformed = pd.read_csv('control_processed.csv', index_col=0).transpose()
    control_transformed['Label'] = 'TD'
    data2 = case_transformed.append(control_transformed).reset_index()
    data2.rename(columns={'index':'Sample Name','Label':'Group'}, inplace=True)
    embed()
    '''
    
     # fill using pcaMetods in R:
    # requires to activate metaboanalyst2 via conda
    # '''
    
    data2 = run_ppca_in_R(data, 'ASD', 'TD')
    
    
        
    # embed()
    # '''
    
    
    # fill from metaboanalyst - remotely:
    '''
    data2 = pd.read_csv('Davis MIND Autism Pre-Vax_M+F_Met+Exp_PROCESSED.csv', index_col=0).transpose().reset_index()
    data2.rename(columns={'index':'Sample Name','Label':'Group'}, inplace=True)
    '''
    
     # fill from metaboanalyst - locally:
    '''
    import os
    data.to_csv('data.csv', index=False)
    
    cmd = 'Rscript.exe C:/Users/Jon/Documents/BITBUCKET/metabotools/R_scripts/ppca_impute.R "%s"'%'data.csv'
    print(cmd)
    os.system(cmd)
    data2 = pd.read_csv('data_processed.csv', index_col=0).transpose().reset_index()
    data2.rename(columns={'index':'Sample Name','Label':'Group'}, inplace=True)
    
    '''
    
    # embed()
    
    # fill using knn
    '''
    data2 = fill_na_knn(data)
    '''
    
    # fill imputed data with red
    data2 = fill_na_red(data, data2)

    
    data2.to_excel('Davis MIND Autism Pre-Vax_M+F_Met+Exp_PPCA.xlsx')
    
def run_ppca_in_R_multiple_groups(data, non_data_cols=2):
    import os
    curdir = os.getcwd()
    os.chdir(r'C:\Users\Jon\Documents\BITBUCKET\metabotools')
    from shutil import move
        
    orig_data = data.copy()
    
    groups = data['Group'].unique()
    data = data.sort_values('Group').reset_index()
    del data['index']
    
    out = {}
    for group in groups:
        data2 = data[data['Group']==group]
        data2.to_csv('temp.csv', index=False)      
        cmd = 'Rscript.exe C:/Users/Jon/Documents/BITBUCKET/metabotools/R_scripts/ppca_impute_pcamethods.R "%s" %i'%('temp.csv', non_data_cols) # +1 because R index start at 1
        print(cmd)
        os.system(cmd)
        out[group] = pd.read_csv('ppca_imputed.csv', index_col=0)
        
    case_transformed = pd.read_csv('case_processed.csv', index_col=0)
    control_transformed = pd.read_csv('control_processed.csv', index_col=0)
    # embed()
    out2 = pd.DataFrame()
    for group in out.keys():
        out2 = out2.append(out[group])
        
    out2 = out2.reset_index()
    del out2['index']
    
    out2 = pd.merge(data[data.columns[:non_data_cols]], out2, left_index=True, right_index=True)
    out2.rename(columns=dict(zip(out2.columns,data.columns)), inplace=True)
    os.chdir(curdir)
    
    out2.index = out2['Sample Name']
    del out2['Sample Name']
    out3 = out2.loc[orig_data['Sample Name'],:] # change order back to original
    out3 = out3.reset_index()
        
    return out3
    
def verify(data, data2, non_data_cols=2):
    # verify that two datasets are the same
    # data.index = data['Sample Name']
    # data2.index = data2['Sample Name']
    res = data[data.columns[non_data_cols:]] - data2[data2.columns[non_data_cols:]]
    # works so long as indices are correct
    return res
    
    
def run_ppca_in_R(data, case, control, non_data_cols=2):
    import os
    # curdir = os.getcwd()
    from shutil import move
    
    # data = data.sort_values('Sample Name').reset_index()
    # del data['index']
    case_data = data[data['Group']==case]
    control_data = data[data['Group']==control]
    
    case_data.to_csv('case.csv', index=False)
    control_data.to_csv('control.csv', index=False)
    
    
    R_script_dir = (os.sep).join(__file__.split(os.sep)[:-1])+os.sep+'R_scripts'
    
    cmd = 'Rscript.exe %s/ppca_impute_pcamethods.R "%s" %i'%(R_script_dir, 'case.csv', non_data_cols) # +1 because R index start at 1
    print(cmd)
    os.system(cmd)
    import time
    time.sleep(2) # moving files too quick causes problems
    move('ppca_imputed.csv','case_processed.csv')
        
        
    cmd = 'Rscript.exe %s/ppca_impute_pcamethods.R "%s" %i'%(R_script_dir, 'control.csv', non_data_cols)
    print(cmd)
    os.system(cmd)
    import time
    time.sleep(2) # moving files too quick causes problems
    
    move('ppca_imputed.csv','control_processed.csv')
    
    case_transformed = pd.read_csv('case_processed.csv', index_col=0)
    case_transformed.index = case_data['Sample Name']
    case_transformed.loc[:,'Group']=case
    control_transformed = pd.read_csv('control_processed.csv', index_col=0)
    control_transformed.index = control_data['Sample Name']
    control_transformed.loc[:,'Group']=control
    data2 = case_transformed.append(control_transformed).reset_index()
    data2 = data2[['Sample Name','Group'] + data2.columns[1:-1].tolist()]
    # del data2['index']
    # data2 = pd.merge(data[data.columns[:non_data_cols]], data2, left_index=True, right_index=True)
    data3 = data2.rename(columns=dict(zip(data2,data.columns)))
    data3.index = data3['Sample Name']
    data4 = data3.loc[data['Sample Name'],:]
    del data4['Sample Name']
    data4 = data4.reset_index()
    # embed()
    
    # os.chdir(curdir)
    return data4
    
    
def drop_too_many_na_cols(data, percent=100):
    num_samples = len(data)
    missing_threshold = num_samples - ((percent/100.)*num_samples)
    print('dropping features with > %i missing values'%missing_threshold)
    cnt = data.count()
    keep = cnt[cnt>missing_threshold]
    data = data[keep.index]
    return data
    
def fill_na_knn(data):
    import numpy as np
    
    ''' doesnt work - experimental '''
    # from sklearn.experimental import enable_iterative_imputer
    # from sklearn.impute import IterativeImputer
    # imp = IterativeImputer(max_iter=10, random_state=0)
    # imp.fit(data)
    # data2 = imp.transform(data)
    
    
    
    ''' knn imputer '''
    from sklearn.impute import KNNImputer
    imputer = KNNImputer(n_neighbors=2, weights="uniform")
    data2 = imputer.fit_transform(data[data.columns[2:]])
    data2 = pd.DataFrame(data2, columns=data.columns[2:])
    data2 = pd.merge(data[data.columns[:2]], data2, left_index=True, right_index=True)
    embed()
    
    return data2
    
def fill_missing(data):
    '''
    This function borrowed from hypertools: https://hypertools.readthedocs.io/en/
    '''
    x = data.values # jon added this
    import numpy as np
    from ppca import PPCA
    # ppca if missing data
    m = PPCA()
    m.fit(data=np.vstack(x)) 
    x_pca = m.transform()

    # if the whole row is missing, return nans
    all_missing = [idx for idx, a in enumerate(np.vstack(x)) if all([type(b)==np.nan for b in a])]
    if len(all_missing)>0:
        for i in all_missing:
            x_pca[i, :] = np.nan
    
    tmp=pd.DataFrame(stats.zscore(data), columns=data.columns)
        
    ''' jon additions '''
    df = data
    df_z = pd.DataFrame(m.data, columns=data.columns) # m.data = z-scores of original values, with new values filled in...
    mean_std={}
    for var in df.columns:
        mean_std[var]=(df[var].mean(), df[var].std())

    def reverse_zscore(pandas_series, mean, std):
        '''Mean and standard deviation should be of original variable before standardization'''
        yis=pandas_series*std+mean
        return yis
        
    embed()
    original_mean, original_std = mean_std[var]
    original_var_series = reverse_zscore(df_z[var], original_mean, original_std)
    
    
    ''' '''
            
            
    embed()
    # get the original lists back
    if len(x)>1:
        
        x_split = np.cumsum([i.shape[0] for i in x][:-1])
        return list(np.split(x_pca, x_split, axis=0))
    else:
        return [x_pca]
    
def fill_na_ppca(data):
    '''
    Uses: https://github.com/el-hult/pyppca
    > pip install git+https://github.com/el-hult/pyppca
    
    #pip install hypertools
    
    
    '''
    import numpy as np
    
    import hypertools as hyp
    
    from ppca import PPCA
    ppca = PPCA()
    
    
    data1 = data[data.columns[2:]]
    
    data2 = fill_missing(data1)
    
    embed()
    
    red = ppca.fit(data1.values)
    red = ppca.fit(np.vstack(data1.values))
    
    # data2= hyp.reduce(data1) # this USES ppca
    
    embed()
    
    data2 = pd.DataFrame(data2, columns=data.columns[2:])
    data2 = pd.merge(data[data.columns[:2]], data2, left_index=True, right_index=True)
    
    embed()
        
    ''' ppca imputer '''
    from pyppca import ppca
    
    only_num = data[data.columns[2:]]
    
    Y = only_num
    d=10
    dia=False
    embed()
    
    C, ss, M, X, Ye = ppca(Y=Y,d=d,dia=dia)
    
    return data2
    
    
def fill_na_col_mean(data):
    '''
    fills nas with the mean of the column, by group
    '''
    print('fillina N/As with mean of column, by group')
    out = {}
    data2 = pd.DataFrame()
    for g in data['Group'].unique():
        tmp = data[data['Group']==g]
        tmp = tmp.fillna(tmp.mean())
        out[g] = tmp
        data2 = data2.append(tmp)
    
    return data2
    
    
def fill_na_red(data, data2): # requires original dataframe and imputed dataframe
    '''
    colors values red if they were na in the original df    
    '''
    
    def apply_color(x):
        # lambda x: np.nan if x < 90 else x
        import math
        def myfunc(x):
            # from IPython import embed
            # embed()
            if type(x)==str:
                return 'background-color: #ffffff'
            elif math.isnan(x):
                return 'background-color: red'
            else:
                return 'background-color: #ffffff'
        
        return data.applymap(lambda val: myfunc(val))
    
    data2 = data2.style.apply(apply_color, axis=None)
    
    
    return data2
    
    #




if __name__=='__main__':
    main()