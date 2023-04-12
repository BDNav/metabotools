import pandas as pd
from IPython import embed


def main():



    # compare()
    
    # embed()


    # read in a "VR" sheet:
    filename_case = 'Suicide_Study_Depress_spearman_VR_SHEET_v1.0.xlsx'
    filename_control = 'Suicide_Study_Control_spearman_VR_SHEET_v1.0.xlsx'
    
     # outputs csv of correlation lines 
    
    build_lines(filename_case)
    
    embed()
    
    r_threshold = 0.5
    FDR_threshold = 1
    # inout = 'OUT'
    inout = ''

    case_sums = load_filter_data(filename_case, r_threshold=r_threshold, FDR_threshold=FDR_threshold, inout=inout)
    control_sums = load_filter_data(filename_control, r_threshold=r_threshold, FDR_threshold=FDR_threshold, inout=inout)

    case_sums.to_csv('case_sums.csv')    
    control_sums.to_csv('control_sums.csv')    
    heights = normalize_heights(case=case_sums, control=control_sums)

    # heights.to_excel('heights.xlsx')
    heights.fillna(0).to_csv('heights.txt', sep=';')
    
    embed()
    
def build_lines(infile):
    data = pd.read_excel(infile,'Unity')
    embed()
    data = data.sort_values('abs(r)', ascending=False)
    out = data[['x1','y1','x2','y2','spearman r']]
    out[:100].to_csv('lines.txt')
    
    

def compare():
    case = pd.read_csv('Nate_sheets_for_comparison/Depress_case_pt5.csv', sep=';')
    control = pd.read_csv('Nate_sheets_for_comparison/Depress_cont_pt5.csv', sep=';')
    
    # heights = pd.read_excel('heights.xlsx')
    
    
    # case = pd.merge(heights, case, left_index=True, right_on='Chem Name')
    # control = pd.merge(heights, control, left_index=True, right_on='Chem Name')
    
    case.to_csv('case_comp.csv')
    control.to_csv('control_comp.csv')
    
    # comp = pd.merge(comp, case, left_on='Chem Name', right_on='Chem Name')
    
    
    embed()


def load_filter_data(filename, r_threshold=0, FDR_threshold=1, Q_threshold=1, inout=''):
    '''
    function to load the data, and filter it based on our chosen thresholds:
    Takes in a VR sheet filename, filters based on thresholds and returns a list of metabolites and their sum(r)
    '''

    data = pd.read_excel(filename,'Unity')

    # we are only interested in the r, fdr, pathway, and a single metabolite (e.g. x1, name1, etc)
    data = data[['spearman r', 'FDR', 'name1', 'Pathway', 'name2']]

    # filter data based on r, FDR and/or in/out of pathway:
    data = data[abs(data['spearman r'])>r_threshold]
    data = data[data['FDR']<FDR_threshold]
    data = data[data['Q value']<Q_threshold]
    if inout!='':
        data = data[data['Pathway']==inout]

    data = data[['name1','spearman r', 'name2']].drop_duplicates()
    # simplify to just a list of metabolites and r scores
    data = data[['name1','spearman r']]
    
    
    # separate pos from neg?
    pos = data[data['spearman r']>=0]
    neg = data[data['spearman r']<0]
    
    pos = pos.groupby('name1').sum().sort_values('spearman r')
    neg = neg.groupby('name1').sum().sort_values('spearman r')
   

    # sum up the r values for each metabolite:
    sums = data.groupby('name1').sum().sort_values('spearman r')
    
    # embed()
    
    return sums
    
def load_filter_data2(data, r_threshold=0, FDR_threshold=1, Q_threshold=1, inout=''):
    '''
    function to load the data, and filter it based on our chosen thresholds:
    Takes in a VR sheet filename, filters based on thresholds and returns a list of metabolites and their sum(r)
    
    Updated to '2' on 9/22/2020 to take in df (not file)
    Also to return pos and negatives (not all)
    
    '''
    print("Calculating heights for r>%0.2f and fdr<%0.2f and Q<%0.2f"%(r_threshold,FDR_threshold,Q_threshold))
    # embed()
    # we are only interested in the r, fdr, pathway, and a single metabolite (e.g. x1, name1, etc)
    data = data[['spearman r', 'FDR', 'name1', 'Pathway', 'name2']]

    # filter data based on r, FDR and/or in/out of pathway:
    data = data[abs(data['spearman r'])>r_threshold]
    data = data[data['FDR']<FDR_threshold]
    data = data[data['FDR']<Q_threshold]
    if inout!='':
        data = data[data['Pathway']==inout]

    data = data[['name1','spearman r', 'name2']].drop_duplicates()
    # simplify to just a list of metabolites and r scores
    data = data[['name1','spearman r']]
    
    # embed()
    # separate pos from neg?
    pos = data[data['spearman r']>=0]
    neg = data[data['spearman r']<0]

    pos_edge_counts = pos.groupby('name1').count()
    neg_edge_counts =  neg.groupby('name1').count()
    pos = pos.groupby('name1').sum().sort_values('spearman r')
    neg = neg.groupby('name1').sum().sort_values('spearman r')
   


    # sum up the r values for each metabolite:
    # sums = data.groupby('name1').sum().sort_values('spearman r')
    
    # embed()
    
    return {'pos':pos, 'neg':neg, 'pos_edges': pos_edge_counts, 'neg_edges': neg_edge_counts}
    
def load_filter_data_lines(filename, r_threshold=.9, FDR_threshold=.1, inout=''):
    '''
    function to load the data, and filter it based on our chosen thresholds:
    Takes in a VR sheet filename, filters based on thresholds and returns a list of metabolite pairs that satisfy thresholds
    '''

    data = pd.read_excel(filename,'Unity')
    # embed()
    # we are only interested in the r, fdr, pathway, and a single metabolite (e.g. x1, name1, etc)
    data = data[['spearman r', 'FDR', 'Chemical Index1', 'Pathway', 'Chemical Index2']]

    # filter data based on r, FDR and/or in/out of pathway:
    data = data[abs(data['spearman r'])>r_threshold]
    data = data[data['FDR']<FDR_threshold]
    if inout!='':
        data = data[data['Pathway']==inout]

    data = data[['Chemical Index1','Chemical Index2','spearman r']].drop_duplicates()
    
    return data


def transform(sums):
    '''
    log2 transforms data, taking into account negative vals and vals < 1
    ''' 
    import numpy as np
    
    pos = sums[sums>0].dropna()
    neg = sums[sums<=0].dropna()
    
    # TEMP: add 1 and -1 to all, to account for < 1 problems
    pos = pos+1
    neg = -1*(neg-1)
    
    pos = np.log2(pos)
    neg = -1*np.log2(neg)
    
    out = pos.append(neg).sort_values(pos.columns[0])
    
    return out


def normalize_heights(case=[], control=[]):
    '''
    function to calculate peak sizes, takes in a list of unique metabolites and sum(r) values:
    log2 transforms all, then normalizes to the maximum peak height
    '''
    
    #log2 transform both sets:
    case_log2 = transform(case)
    control_log2 = transform(control)
    
    
    # get absolute max value:
    case_max = abs(case_log2).max()[0]
    control_max = abs(control_log2).max()[0]
    total_max = max(case_max, control_max)
       
    case = case_log2/total_max
    case.rename(columns={'spearman r':"case"}, inplace=True)
    control = control_log2/total_max
    control.rename(columns={'spearman r':"control"}, inplace=True)
    out = pd.concat([case, control],axis=1) # note that NANs come from those that don't pass thresholds in one cond or another

    return out

if __name__=='__main__':
    main()



