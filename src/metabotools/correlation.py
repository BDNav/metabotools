

import pandas as pd
from IPython import embed
from scipy import stats
import numpy as np

def get_zscores(data, md_col, case, control,  log_transform=True, filter_missing=True):
    
    data2=pd.DataFrame()
    for c in data.columns[md_col:]:
        # try:
        # print c
        data2[c]=pd.to_numeric(data[c])
        # except:
            # pass
    # data2 = data2.transpose().dropna().transpose()
    # embed()
    # Log transform data
    if log_transform:
        log2 = np.log2(data2)
        log2['Group']=data['Group']
        data = log2
        print("NOTE we are Log2 ing this data, is that what you want?")
        # embed()
    else:
        data2['Group'] = data['Group']
        data = data2
    # embed()
    # split case and controls
    control = data[data['Group']==control][data.columns[:-1]]
    case = data[data['Group']==case][data.columns[:-1]]

    # calculate z-scores
    control_z = pd.DataFrame(stats.zscore(control, ddof=1), index=control.index, columns=control.columns)
    control_mean = control.mean()
    control_std = control.std()
    
    # embed()
    
    case_z = (case-control_mean)/control_std
    # embed()
    # only return cols with data:
    if filter_missing==True:
        show=[]
        for c in data.columns[:-1]:
            if str(data[c].mean())!='nan':
                show.append(c)
    else:
        show = case_z.columns
    # embed()
            
    # control_z[show].to_csv('GWI_Control_zscore.csv')
    # case_z[show].to_csv('GWI_GWS_zscore.csv')

    
    return {'case':case_z[show], 'control':control_z[show]}