from IPython import embed
    
import multiprocessing as mp
# from IPython import embed
import pandas as pd

def main():
    # TESTING:
    
    test_run_bootstrap()
    import os
    os.chdir(r'D:\Dropbox (Personal)\naviauxlab_informatics\ASD_TD\Sai_Correlations_080621')
    from glob import glob
    bootstrap_files = glob('resultsFeb11/*.xlsx')
    analyze_bootstraps(bootstrap_files)
    
def organize_raw_data(dir, name, method, type):

    files = glob(dir + '%s*%s*.xlsx'%(type,method))
    out = []
    dupes_out = {}
    all_samples = {}
    for file in files:
        print(file)
        subsample_size = int(file.split('_')[3])
        data = pd.read_excel(file,'Significant samples')
        
        for c in data.columns[6:]:
            sig_corrs = data[c].sum()
            if c!='full_set':
                out.append({'subsample_size':subsample_size, 'sig_corrs':sig_corrs,'bootstrap_iteration':int(c.split('_')[-1])})
            else:
                out.append({'subsample_size':c, 'sig_corrs':sig_corrs,'bootstrap_iteration':1})
            
        samples = pd.read_excel(file,type+'S')
        samples = samples[samples.columns[2:]]
        # dupes = samples[samples.columns[1:]]
        tmp = samples-1
        dupes = tmp.sum()
        # embed()
        
        dupes_out[subsample_size] = dupes
        all_samples[subsample_size] = samples
        
    dupes_out = pd.DataFrame(dupes_out)
            
    out = pd.DataFrame(out)
    res = pd.pivot_table(out, index='bootstrap_iteration', columns='subsample_size', values='sig_corrs')

    res.loc['Mean',:] = res.mean()
    res.loc['Median',:] = res.median()
    res.loc['Min',:] = res.min()
    res.loc['Max',:] = res.max()
    res.loc['SD',:] = res.std()
    res.loc['SEM',:] = res.sem()
    
    
    res.to_excel(f'{name}_{type}_all_sig_pairs_raw_data.xlsx')
    dupes_out.to_excel(f'{name}_{type}_duplicates.xlsx')
    
    writer = pd.ExcelWriter(f'{name}_{type}_sample_counts.xlsx')
    for s in all_samples.keys():
        all_samples[s].to_excel(writer,sheet_name=str(s))
    
    writer.close()
    
def test_run_bootstrap():
    
    DEBUG=True
    
    import os
    os.chdir(r'D:\Dropbox (Personal)\naviauxlab_informatics\ASD_TD\Sai_Correlations_080621')
    
    data_file = 'Davis_ASD_TD_Pre_endogenous_toxins_75percent_merged_New_20210726_PPCA.csv' # NOTE: deleted index + extra (non met) columns and saved xlsx file to csv for corr analysis
    
    data = data = pd.read_csv(data_file, index_col=0)
    data = data[data.columns[:-19]] # remove toxins so we are only looking at mets    

    if DEBUG==True:
        data = data[data.columns[:25]]

    method = 'spearman'
    tableS1_file = '5-ASD_pre_vs_TD_post_ASD_Pre_TableS1.xlsx' # NOTE THIS MIGHT BE THE WRONG TS1 - double check
    db_version = '7.156'
    db_file = r'MasterDatabase_v%s.xlsx'%db_version
    exp = '2_ASD_Pre'
    control = '1_TD_Pre'
    
    OUTDIR='testDec62022' # where files will be saved

    size = 5
    max_bootstraps = 10
    use_zscore = True
    log_transform = True
    leavein=True
    # (size, casecon, OUTDIR, data, data_file, tableS1_file, exp, control, db_file, method, max_bootstraps, use_zscore, log_transform, Q_threshold=0.05, p_threshold=None, r_threshold=None, leavein=True)
    run_bootstraps(size, 'CONTROLS', OUTDIR, data, data_file, tableS1_file, exp, control, db_file, method, max_bootstraps, use_zscore, log_transform, leavein)
    
    
def analyze_bootstraps(name, method, input_files, data_file, tableS1_file, db_file, transformation):
    '''
    Takes an input of bootstrap files, must have "CASE" and "CONTROL" in name
    More structure reqs:
    
    
    '''
    import datetime as dt
    import pandas as pd
    from scipy import stats
        
    tables1 = pd.read_excel(tableS1_file, 'TableS1')
    knnClusters = pd.read_excel(tableS1_file, 'kNN (10) Ranking by MDA')
    
    
    pwy = pd.read_excel(db_file, 'Pathway', index_col=0)
    chem = pd.read_excel(db_file, 'Chemical')

    alias = pd.read_excel(db_file, 'Alias', index_col=0)
    mrm = pd.read_excel(db_file, 'MRM', index_col=0)

    mtdb = pd.merge(pwy, chem, left_on='Pathway Index', right_on='Pathway Index')
    pwy_count = mtdb.groupby('Pathway Name').count()

    orig_data = pd.read_csv(data_file, index_col=0)

    mets = pd.read_csv(data_file, index_col=0)
    mets = pd.DataFrame(mets.columns[1:])
    mets = pd.merge(alias, mets, left_on='MRM Name', right_on=0)
    measured_mets = mets[['MRM Index','MRM Name']]

    measured_mets = pd.merge(measured_mets, mrm, left_on='MRM Index', right_on='MRM Index')

    measured_chemicals = measured_mets['Chemical Index'].unique()
    measured_pathways = mtdb[mtdb['Chemical Index'].isin(measured_chemicals)]

    met_per_pwy = measured_pathways.groupby('Pathway Name').count()['Chemical Index']
    
    
    mets = pd.read_csv(data_file, index_col=0)
    mets = pd.DataFrame(mets.columns[1:])
    mets = pd.merge(alias, mets, left_on='MRM Name', right_on=0)
    measured_mets = mets[['MRM Index','MRM Name']]
    
    
    metrics = pd.read_csv(data_file, index_col=0)
    metrics2 = pd.read_excel(input_files[0],'Header', index_col=0) # steal header info from one of the bootstraps
    
    # FINDING which label is for case and which is for control
    if 'CASES' in input_files[0]:
        case = metrics2.loc['Cohort Description--Diagnosis (Case, Control?)'][0]
        control = list(set(metrics['Group'].unique()) - set([case]))[0]
    else:
        control = metrics2.loc['Cohort Description--Diagnosis (Case, Control?)'][0]
        case = list(set(metrics['Group'].unique()) - set([control]))[0]
    case_size = len(metrics[metrics['Group']==case])
    control_size = len(metrics[metrics['Group']==control])
    
    targeted_mets_count = len(orig_data.columns)-1
    my_mets = orig_data.transpose().dropna().transpose()
    measured_mets_count = len(my_mets.columns)

    
    header = {'Info':{
    'Study name':name,
    'Input file':data_file,
    'TableS1 file':tableS1_file,
    'MTDB':db_file,
    'method':method,
    'date': dt.datetime.today().strftime("%m/%d/%Y"),
    'Control group name for Z score calculation': control,
    'Case group name for Z score calculation': case,
    'Sample size of control group, N =': control_size,
    'Sample size of case group, N=':case_size,
    'Number of targeted metabolites':targeted_mets_count,
    'Number of measured metabolites':measured_mets_count,
    'Total possible combinations, K = ':(measured_mets_count**2-measured_mets_count)/2,
    'transformation': transformation
    }}
    
    header = pd.DataFrame(header)
    
    # embed()
    # get mets per pathways (both canonical and knn):  build dict of pathways: list of mets
    pwy_mets = {}
    for up in tables1['K-Means Cluster (of 10)'].unique():
        mets = tables1[tables1['K-Means Cluster (of 10)']==up]['MRM Name'].tolist()
        pwy_mets['kNN cluster ' + str(up)]=mets
    for up in tables1['Pathway Name'].unique():
        mets = tables1[tables1['Pathway Name']==up]['MRM Name'].tolist()
        pwy_mets[up]=mets
    # embed()
    # end mapping metes

    out = []
    pwys_out = []
    nodes_out = []
    pwy_nodes_out = []
    all_rstats = []
    all_pwy_rstats = []

    def calc_stats(data, exp, size):
        res = []
        mean = data.sum().mean()
        median = data.sum().median()
        std = data.sum().std()
        min = data.sum().min()
        max = data.sum().max()
        lci, hci = stats.norm.interval(alpha=0.95, loc=data.sum().mean(), scale=data.sum().sem())
        res = {
            'type':exp,
            'size':size,
            'mean':mean,
            'median':median,
            'std':std,
            'min':min,
            'max':max,
            'L_ci':lci,
            'H_ci':hci
        }
        
        return res
        
        
    def calc_node_stats(data, exp, size):
        node_counts = {}
        for c in data.columns[5:-1]:
            # print(c)
            tmp = data[['from','to',c]].dropna()
            if len(tmp)>0:
                nodes = len(set(tmp['from'].tolist()+tmp['to'].tolist()))
            else:
                nodes = 0
            node_counts[c] = nodes
            
        node_counts = pd.Series(node_counts)
        
        res = []
        mean = node_counts.mean()
        median = node_counts.median()
        std = node_counts.std()
        min = node_counts.min()
        max = node_counts.max()
        lci, hci = stats.norm.interval(alpha=0.95, loc=node_counts.mean(), scale=node_counts.sem())
        
        res = {
            'type':exp,
            'size':size,
            'Nodes mean':mean,
            'Nodes median':median,
            'Nodes std':std,
            'Nodes min':min,
            'Nodes max':max,
            'Nodes L_ci':lci,
            'Nodes H_ci':hci
        }
        return res
            
        
    def get_all_stats(data, sign, is_max=False):
        # unique_pathways = list(set(data['Pathway1'].unique().tolist() + data['Pathway2'].unique().tolist()))
        unique_pathways = pwy_mets.keys()
        out = []
        nodes_out = []
        pwy_nodes_out = []
        pwys_out = []
        for pwy in unique_pathways: # pathways defined by table s1
            mets = pwy_mets[pwy]
            data2 = data[data['from'].isin(mets)]
            data3 = data[data['to'].isin(mets)]
            umets = data2['from'].unique().tolist()
            umets += data3['to'].unique().tolist()
            
            # data2 = data[data['Pathway1']==pwy]
            # umets = data2['from'].unique().tolist()
            # data3 = data[data['Pathway2']==pwy]
            # umets += data3['to'].unique().tolist()
            
            data4 = pd.concat([data2,data3])
            node_res = calc_node_stats(data4, exp, size)
            node_res['sign']=sign
            node_res['pathway']=pwy
            pwy_nodes_out.append(node_res)
            data4 = data4[data4.columns[5:-1]]
           
            pwy_res = calc_stats(data4, exp, size)
            pwy_res['pathway']=pwy
            pwy_res['sign']=sign
            pwys_out.append(pwy_res)

        for ptype in data['Pathway'].unique(): # in out or both
            data2 = data[data['Pathway']==ptype]
            node_res = calc_node_stats(data2, exp, size)
            node_res['ptype']=ptype
            node_res['sign']=sign
            nodes_out.append(node_res)
            data2 = data2[data2.columns[5:-1]]
            
            res = calc_stats(data2, exp, size)
            res['ptype']=ptype
            res['sign']=sign
            out.append(res)
        
        # in and out:
        node_res = calc_node_stats(data, exp, size)
        node_res['ptype']='ALL'
        node_res['sign']=sign
        nodes_out.append(node_res)
        data2 = data[data.columns[5:-1]]
        res = calc_stats(data2, exp, size)
        res['ptype']='ALL'
        res['sign']=sign
        out.append(res)
        
        
        return out, pwys_out, nodes_out, pwy_nodes_out
        
    def get_rstats(data, exp, size):
        all_res = []
        data = data[data['bootstrap']!='full_set']
        for type in data['type'].unique(): # all, positive, negative
            data2 = data[data['type']==type]
            
            for ptype in data2['ptype'].unique():
                data3 = data2[data2['ptype']==ptype]
                lci, hci = stats.norm.interval(alpha=0.95, loc=data3['mean r'].mean(), scale=data3['mean r'].sem())
                res = {
                    'type':exp,
                    'size':size,
                    'r-value mean':data3['mean r'].mean(), # note: mean of a means...
                    'r-value median':data3['median r'].mean(),  # mean medians
                    'r-value std':data3['std r'].mean(),
                    'r-value min':data3['min r'].min(),
                    'r-value max':data3['max r'].max(),
                    'r-value L_ci':lci,
                    'r-value H_ci':hci,
                    'sign': type, # positive, negative, all
                    'ptype':ptype
                }
                all_res.append(res)
        return all_res
        
    def get_pwy_rstats(data, exp, size):
        all_res = []
        for pwy in data['pathway'].unique():
            data2 = data[data['pathway']==pwy]
            
            lci, hci = stats.norm.interval(alpha=0.95, loc=data2['mean r'].mean(), scale=data2['mean r'].sem())
            res = {
                'type':exp,
                'size':size,
                'r-value mean':data2['mean r'].mean(), # note: mean of a means...
                'r-value median':data2['median r'].mean(),  # mean medians
                'r-value std':data2['std r'].mean(),
                'r-value min':data2['min r'].min(),
                'r-value max':data2['max r'].max(),
                'r-value L_ci':lci,
                'r-value H_ci':hci,
                'pathway':pwy
            }
            all_res.append(res)
        return all_res
        
    ''' CALCULATING FULL STATS '''
    file_stats = []
    # get max size for full stat calcs
    for file in input_files:

        size = int(file.split('_')[3]) # e.g. 10, 12, 14, etc
        
        exp = file.split('/')[-1].split('_')[0] # CASES or CONTROLS
        
        file_stats.append({'size':size,'type':exp,'file':file})
    size = 'FULL SET ()'
    file_stats = pd.DataFrame(file_stats)
    
    max_case = file_stats[file_stats['type']=='CASES'].sort_values('size')['file'].tolist()[-1]
    max_control = file_stats[file_stats['type']=='CONTROLS'].sort_values('size')['file'].tolist()[-1]
    
    exp = 'CASE'
    size = 'FULL SET'
    data = pd.read_excel(max_case, 'Significant samples', index_col=0)
    
    data = data[data.columns[:5].tolist()+['full_set',data.columns[6]]]
    pos = pd.read_excel(max_case, 'POS Significant samples', index_col=0)
    pos = pos[pos.columns[:5].tolist()+['full_set',pos.columns[6]]]
    neg = pd.read_excel(max_case, 'NEG Significant samples', index_col=0)
    neg = neg[neg.columns[:5].tolist()+['full_set',neg.columns[6]]]
    
    
    res, pwys, nodes, pwy_nodes = get_all_stats(data, 'all')
    out+=res
    pwys_out+=pwys
    nodes_out+=nodes
    pwy_nodes_out+= pwy_nodes
    
    res, pwys, nodes, pwy_nodes = get_all_stats(pos, 'positive')
    out+=res
    pwys_out+=pwys
    nodes_out+=nodes
    pwy_nodes_out+= pwy_nodes
    
    res, pwys, nodes, pwy_nodes = get_all_stats(neg, 'negative')
    out+=res
    pwys_out+=pwys
    nodes_out+=nodes
    pwy_nodes_out+= pwy_nodes
    
    exp = 'CONTROL'
    size = 'FULL SET'
    data = pd.read_excel(max_control, 'Significant samples', index_col=0)
    
    data = data[data.columns[:5].tolist()+['full_set',data.columns[6]]]
    pos = pd.read_excel(max_control, 'POS Significant samples', index_col=0)
    pos = pos[pos.columns[:5].tolist()+['full_set',pos.columns[6]]]
    neg = pd.read_excel(max_control, 'NEG Significant samples', index_col=0)
    neg = neg[neg.columns[:5].tolist()+['full_set',neg.columns[6]]]
    
    
    res, pwys, nodes, pwy_nodes = get_all_stats(data, 'all')
    out+=res
    pwys_out+=pwys
    nodes_out+=nodes
    pwy_nodes_out+= pwy_nodes
    
    res, pwys, nodes, pwy_nodes = get_all_stats(pos, 'positive')
    out+=res
    pwys_out+=pwys
    nodes_out+=nodes
    pwy_nodes_out+= pwy_nodes
    
    res, pwys, nodes, pwy_nodes = get_all_stats(neg, 'negative')
    out+=res
    pwys_out+=pwys
    nodes_out+=nodes
    pwy_nodes_out+= pwy_nodes
    
    ''' END CALCULATE FULL STATS '''
        

    for file in input_files:

        size = int(file.split('_')[3]) # e.g. 10, 12, 14, etc
        
        exp = file.split('/')[-1].split('_')[0] # CASES or CONTROLS
        
        print(size,exp)
        
        # N stats - node based
        data = pd.read_excel(file, 'Significant samples', index_col=0)
        pos = pd.read_excel(file, 'POS Significant samples', index_col=0)
        neg = pd.read_excel(file, 'NEG Significant samples', index_col=0)
        
        res, pwys, nodes, pwy_nodes = get_all_stats(data, 'all')
        out+=res
        pwys_out+=pwys
        nodes_out+=nodes
        pwy_nodes_out+= pwy_nodes
        
        res, pwys, nodes, pwy_nodes = get_all_stats(pos, 'positive')
        out+=res
        pwys_out+=pwys
        nodes_out+=nodes
        pwy_nodes_out+= pwy_nodes
        
        res, pwys, nodes, pwy_nodes = get_all_stats(neg, 'negative')
        out+=res
        pwys_out+=pwys
        nodes_out+=nodes
        pwy_nodes_out+= pwy_nodes
        
        all_stats = pd.read_excel(file, 'Stats', index_col=0)

        # TMP added because ptype was 'all' intead of 'ALL' 4/25/2022
        toALL = all_stats[all_stats['ptype']=='all'].index
        all_stats.loc[toALL,'ptype']='ALL'
        # end TMP
        
        rstats = get_rstats(all_stats, exp, size)
        all_rstats+=rstats
        
        pwy_stats = pd.read_excel(file, 'ALL PWY STATS', index_col=0)
        pwy_rstats = get_pwy_rstats(pwy_stats, exp, size)
        all_pwy_rstats+=pwy_rstats
        
    out = pd.DataFrame(out).fillna(0)
    out = out.sort_values(['type','size'], ascending=True)

    nodes_out = pd.DataFrame(nodes_out).fillna(0)
    nodes_out = nodes_out.sort_values(['type','size'], ascending=True)

    all_rstats = pd.DataFrame(all_rstats).fillna(0)
    all_rstats = all_rstats.sort_values(['type','size'], ascending=True)

    out2 = pd.merge(out, nodes_out, how='left')

    out3 = pd.merge(out2, all_rstats, how='left')

    writer = pd.ExcelWriter('%s_bootstrap_results_%s.xlsx'%(name,method))
    header.to_excel(writer,'Header')
    for sign in out3['sign'].unique():
        out3[out3['sign']==sign].to_excel(writer,sign)
    
    
    
    writer.save()



    pwys_out = pd.DataFrame(pwys_out).fillna(0)
    # pwys_out['type']=pwys_out['type'].replace('SAMPLES','CASES')
    pwys_out = pwys_out.sort_values(['type','size'], ascending=True)

    pwy_nodes_out = pd.DataFrame(pwy_nodes_out).fillna(0)
    pwy_nodes_out = pwy_nodes_out.sort_values(['type','size'], ascending=True)

    all_pwy_rstats = pd.DataFrame(all_pwy_rstats).fillna(0)
    all_pwy_rstats['sign'] = 'all'

    pwy_stats_out = pd.merge(pwys_out, pwy_nodes_out)
    pwy_stats_out = pd.merge(pwy_stats_out, all_pwy_rstats, how='left')


    samples = pwy_stats_out[pwy_stats_out['type']=='CASES']
    controls = pwy_stats_out[pwy_stats_out['type']=='CONTROLS']

    max_samples = samples['size'].max()
    max_controls = controls['size'].max()

    sample_pwys_out = samples[samples['size']==max_samples]
    control_pwys_out = controls[controls['size']==max_controls]



    writer = pd.ExcelWriter('%s_bootstrap_results_by_pathway_%s.xlsx'%(name,method))
    header.to_excel(writer,'Header')
    knnClusters.to_excel(writer, 'knn clusters')
    control_pwys_out[control_pwys_out['sign']=='all'].to_excel(writer, 'Control Pathway Summary ALL')
    sample_pwys_out[sample_pwys_out['sign']=='all'].to_excel(writer, 'Cases Pathway Summary ALL')
    
    control_pwys_out[control_pwys_out['sign']=='positive'].drop(columns=['r-value mean','r-value median','r-value std','r-value min','r-value max','r-value L_ci','r-value H_ci']).to_excel(writer, 'Control Pathway Summary POS')
    sample_pwys_out[sample_pwys_out['sign']=='positive'].drop(columns=['r-value mean','r-value median','r-value std','r-value min','r-value max','r-value L_ci','r-value H_ci']).to_excel(writer, 'Cases Pathway Summary POS')

    control_pwys_out[control_pwys_out['sign']=='negative'].drop(columns=['r-value mean','r-value median','r-value std','r-value min','r-value max','r-value L_ci','r-value H_ci']).to_excel(writer, 'Control Pathway Summary NEG')
    sample_pwys_out[sample_pwys_out['sign']=='negative'].drop(columns=['r-value mean','r-value median','r-value std','r-value min','r-value max','r-value L_ci','r-value H_ci']).to_excel(writer, 'Cases Pathway Summary NEG')
        
    writer.save()


def run_parallel_bootstraps(OUTDIR, data_file, tableS1_file, exp, control, db_file, method, start, step_size, max, use_zscore=True, log_transform=True, leavein=False):
    
    from functools import partial
    from itertools import repeat
    from multiprocessing import Pool, freeze_support

    data = data = pd.read_csv(data_file, index_col=0)
    samples = data[data['Group']==exp].index.tolist()
    controls = data[data['Group']==control].index.tolist()

    # for casecon in ['CONTROLS']:
    for casecon in ['CASES','CONTROLS']:
        if casecon=='CONTROLS':
            sizes = range(start,len(controls)+1,step_size) 
        elif casecon=='CASES':
            sizes = range(start,len(samples)+1,step_size)
        else:
            print('error, invalid type specified')
        pool = Pool()
        pool.starmap(run_bootstraps, zip(sizes, repeat(casecon), repeat(OUTDIR), repeat(data), repeat(data_file), repeat(tableS1_file), repeat(exp), repeat(control), repeat(db_file), repeat(method), repeat(max), repeat(use_zscore), repeat(log_transform), repeat(leavein)))
    
    
    
def run_bootstraps(size, casecon, OUTDIR, data, data_file, tableS1_file, exp, control, db_file, method, max_bootstraps, use_zscore, log_transform, leavein, Q_threshold=0.05, p_threshold=None, r_threshold=None):
    '''
    type = 'CASES' or 'CONTROLS'
    '''
    
    # from metabotools.tableOne import *
    # from IPython import embed
    from datetime import date
    from metabotools.correlation_analysis import correlation_analysis
    import random

    today = date.today()

    tables1 = pd.read_excel(tableS1_file, 'TableS1')

    num_metabolites = len(data.columns)-1
    samples = data[data['Group']==exp].index.tolist()
    controls = data[data['Group']==control].index.tolist()

    # calculating total possible combinations (factorial see 10_bootstraps_11192021_rkn2.xlsx
    from math import factorial
    n_cases = len(samples)
    n_controls = len(controls)
    if n_cases>n_controls:
        tot_possible = factorial(n_cases) / (factorial(n_controls) * factorial(n_cases-n_controls))
    else:
        tot_possible = factorial(n_cases) / (factorial(n_controls) * factorial(n_controls-n_cases))


    header = {'value':{
        'Data csv file used': data_file, 
        'Date of Analysis': today.strftime("%B %d, %Y"),
        'Cohort Description--Diagnosis (Case, Control?)': exp,
        'Total number of samples':len(samples),
        'Sample number used for each bootstrap':len(controls),
        'Total possible combinations': format(tot_possible, "5.7E"),
        'Correlation method':method,
        'Group 1 AUC types for correlation (metabolites or toxins?':'metabolites',
        'Group 2 AUC types for correlation (metabolites or toxins?':'metabolites',
        'Number of Group 1 parameters to be correlated':num_metabolites,
        'Number of Group 2 parameters to be correlated':num_metabolites,
        'Total possible unique correlation pairs':(num_metabolites**2-num_metabolites)/2,
        'q threshold for filtering': '<= %0.2f'%Q_threshold,
        'r threshold for filtering': r_threshold,
        'p threshold for filtering': p_threshold,
        'In/Out/All correlations': 'All',
        'log_tranform': log_transform,
        'zscores': use_zscore
    }}

    
    controls = data[data['Group']==control].index.tolist()
    cases = data[data['Group']==exp].index.tolist()

    if leavein==False: # original bootstrapping method - pull samples out, no dupes
        from itertools import combinations
        use_combos = False
        
        if casecon=='CONTROLS':
            if size>n_controls-3:
                combos = [list(c) for c in combinations(controls, size)]
                use_combos=True
                random.shuffle(combos)
                combos.insert(0,controls)
            
        elif casecon=='CASES':
            if size>n_cases-3:
                combos = [list(c) for c in combinations(cases, size)]
                use_combos=True
                random.shuffle(combos)
                combos.insert(0,cases)
    else:  # new bootstrapping method - leave samples in, no need to iterate possible combos
        use_combos = False
    
    orig_max_bootstraps = max_bootstraps

    out = pd.DataFrame()
    pos_out = pd.DataFrame()
    neg_out = pd.DataFrame()
    bootstrap_samples = pd.DataFrame()
    stats = []
    all_pwy_stats = []

    if use_combos==True and max_bootstraps>len(combos):
        max_bootstraps = len(combos)-1
    print('running %i bootstraps for %s with size %i'%(max_bootstraps, casecon, size))
    
    for i in range(0,max_bootstraps+1):
        
        if leavein==False: # original bootstrapping method - pull samples out, no dupes
            if casecon=='CONTROLS':
                subset = random.sample(controls, size)
                if i == 0:
                    subset = controls
            elif casecon=='CASES':
                subset = random.sample(cases, size)
                if i == 0:
                    subset=cases
        else: # bootstrapping with leaveing - e.g. random sampling with replacement
            import numpy as np
            if casecon=='CONTROLS':
                subset = list(np.random.choice(controls, size)) # bootstrapping with leaveing - e.g. random sampling with replacement
                if i == 0:
                    subset = controls
            elif casecon=='CASES':
                subset = list(np.random.choice(cases, size))
                if i == 0:
                    subset=cases
            
        # if there are fewer combos than the boostraps:
        if use_combos==True:
            subset=combos[i]

        print('bootstrap %i for %i samples'%(i, len(subset)))
        
        bootstrap_samples = pd.concat([bootstrap_samples, (pd.DataFrame(['boostrap_%i'%i]*len(subset), subset))])
        

        if casecon=='CONTROLS':
            data2 = data.loc[subset+samples,:] # subset of controls
        elif casecon=='CASES':
            data2 = data.loc[controls+subset,:] # subset of samples
        
        if len(subset)!=len(set(subset)): # need to generate new indices so that sample indexes are unique
            print('dupes')
            data2['pos'] = range(0,len(data2))
            data2['new_index'] = data2['pos'].astype(str) + '_' + data2.index
            data2.index = data2['new_index']
            data2 = data2.drop(columns=['pos','new_index'])
            data2.index.name = 'Sample Name'
            
        
        corrs = correlation_analysis.bootstrap(data2, control, exp, db_file, method=method, use_zscore=use_zscore, log_transform=log_transform) # this will now be TableS2
        
        # if len(subset)!=len(set(subset)):
        #     embed()
        
        if casecon=='CONTROLS':
            corr_results = corrs['%s_control'%method]
        elif casecon=='CASES':
            corr_results = corrs['%s_case'%method]
            
        # FILTERING:
        
        corr_results = corr_results[corr_results['Q value'] < Q_threshold]
        
        def get_stats(input):
            stats = {
            'median p':input['p'].median(),
            'median Q':input['Q value'].median(),
            'median r':input['r'].median(),
            'max p':input['p'].max(),
            'max Q':input['Q value'].max(),
            'max r':input['r'].max(),
            'min p':input['p'].min(),
            'min Q':input['Q value'].min(),
            'min r':input['r'].min(),
            'mean r': input['r'].mean(),
            'std r': input['r'].std()
            }
            return stats

        if i==0:
            cur_bootstrap = 'full_set'
        else:   
            cur_bootstrap = 'boostrap_%i'%i
            
        # collect all corrs:
        in_corr = corr_results[corr_results['Pathway']=='IN']
        out_corr = corr_results[corr_results['Pathway']=='OUT']
        neg_corr = corr_results[corr_results['r']<0]
        pos_corr = corr_results[corr_results['r']>=0]

        in_pos_corr = in_corr[in_corr['r']>=0]
        in_neg_corr = in_corr[in_corr['r']<0]
        out_pos_corr = out_corr[out_corr['r']>=0]
        out_neg_corr = out_corr[out_corr['r']<0]

        # collect all stats:
        all_stats = get_stats(corr_results)
        all_stats['type'] = 'all'
        all_stats['ptype'] = 'ALL'
        all_stats['bootstrap'] = cur_bootstrap
        
        in_stats = get_stats(in_corr)
        in_stats['type'] = 'all'
        in_stats['ptype'] = 'IN'
        in_stats['bootstrap'] = cur_bootstrap

        out_stats = get_stats(out_corr)
        out_stats['type'] = 'all'
        out_stats['ptype'] = 'OUT'
        out_stats['bootstrap'] = cur_bootstrap
        
        pos_stats = get_stats(pos_corr)
        pos_stats['type'] = 'positive'
        pos_stats['ptype'] = 'ALL'
        pos_stats['bootstrap'] = cur_bootstrap

        neg_stats = get_stats(neg_corr)
        neg_stats['type'] = 'negative'
        neg_stats['ptype'] = 'ALL'
        neg_stats['bootstrap'] = cur_bootstrap

        in_pos_stats = get_stats(in_pos_corr)
        in_pos_stats['type'] = 'positive'
        in_pos_stats['ptype'] = 'IN'
        in_pos_stats['bootstrap'] = cur_bootstrap

        in_neg_stats = get_stats(in_neg_corr)
        in_neg_stats['type'] = 'negative'
        in_neg_stats['ptype'] = 'IN'
        in_neg_stats['bootstrap'] = cur_bootstrap

        out_pos_stats = get_stats(out_pos_corr)
        out_pos_stats['type'] = 'positive'
        out_pos_stats['ptype'] = 'OUT'
        out_pos_stats['bootstrap'] = cur_bootstrap

        out_neg_stats = get_stats(out_neg_corr)
        out_neg_stats['type'] = 'negative'
        out_neg_stats['ptype'] = 'OUT'
        out_neg_stats['bootstrap'] = cur_bootstrap

        stats.append(all_stats)
        stats.append(in_stats)
        stats.append(out_stats)
        stats.append(pos_stats)
        stats.append(neg_stats)
        stats.append(in_pos_stats)
        stats.append(in_neg_stats)
        stats.append(out_pos_stats)
        stats.append(out_neg_stats)
        
        
        # getting canonical pathway stats:
        for unique_pwy in tables1['Pathway Name'].unique():
            mets_in_pwy = tables1[tables1['Pathway Name']==unique_pwy]['MRM Name'].tolist()
            pwy_cors = corr_results[corr_results['from'].isin(mets_in_pwy)]            
            pwy_cors = pd.concat([pwy_cors,corr_results[corr_results['to'].isin(mets_in_pwy)]])       
            pwy_stats = get_stats(pwy_cors)
            pwy_stats['pathway']=unique_pwy
            pwy_stats['bootstrap']='boostrap_%i'%i
            all_pwy_stats.append(pwy_stats)
        
        # getting knn pathway stats:
        for unique_pwy in tables1['kNN Cluster (of 10)'].unique():
            mets_in_pwy = tables1[tables1['kNN Cluster (of 10)']==unique_pwy]['MRM Name'].tolist()
            pwy_cors = corr_results[corr_results['from'].isin(mets_in_pwy)]            
            pwy_cors = pd.concat([pwy_cors,corr_results[corr_results['to'].isin(mets_in_pwy)]])       
            pwy_stats = get_stats(pwy_cors)
            pwy_stats['pathway']='kNN cluster ' + str(unique_pwy)
            pwy_stats['bootstrap']='boostrap_%i'%i
            all_pwy_stats.append(pwy_stats)
            
        
        pos_corr_results = corr_results[corr_results['r'] >=0]
        neg_corr_results = corr_results[corr_results['r'] < 0]
        
        if i==0:
            if len(corr_results)>0:
                corr_results.loc[:,'bootstrap'] = 'full_set'
            if len(pos_corr_results)>0:
                pos_corr_results.loc[:,'bootstrap'] = 'full_set'
            if len(neg_corr_results)>0:
                neg_corr_results.loc[:,'bootstrap'] = 'full_set'
        else:
            if len(corr_results)>0:
                corr_results.loc[:,'bootstrap'] = 'boostrap_%i'%i
            if len(pos_corr_results)>0:
                pos_corr_results.loc[:,'bootstrap'] = 'boostrap_%i'%i
            if len(neg_corr_results)>0:
                neg_corr_results.loc[:,'bootstrap'] = 'boostrap_%i'%i
            
        out = pd.concat([out, corr_results])
        pos_out = pd.concat([pos_out, pos_corr_results])
        neg_out = pd.concat([neg_out, neg_corr_results])
                    
        
    out['pres'] = 1
    pos_out['pres'] = 1
    neg_out['pres'] = 1
    matrix=pd.pivot_table(out, index=['from','to','Pathway1', 'Pathway2', 'Pathway'], columns='bootstrap', values = 'pres')
    matrix = matrix.reset_index()
    matrix_pos=pd.pivot_table(pos_out, index=['from','to','Pathway1', 'Pathway2', 'Pathway'], columns='bootstrap', values = 'pres')
    matrix_pos = matrix_pos.reset_index()
    matrix_neg=pd.pivot_table(neg_out, index=['from','to','Pathway1', 'Pathway2', 'Pathway'], columns='bootstrap', values = 'pres')
    matrix_neg = matrix_neg.reset_index()

    bootstrap_samples_out = bootstrap_samples.reset_index()
    bootstrap_samples_out['pres']=1
    bootstrap_samples_out = bootstrap_samples_out.rename(columns={0:'bootstrap'})
    samples_out = pd.pivot_table(bootstrap_samples_out, index='index', columns='bootstrap', values='pres', aggfunc=sum)

    all_pwy_stats = pd.DataFrame(all_pwy_stats)


    out_writer = pd.ExcelWriter(OUTDIR+'/%s_%i_bootstraps_%i_samples_%s.xlsx'%(casecon,orig_max_bootstraps,size,method))
    pd.DataFrame(header).to_excel(out_writer, 'Header')
    pd.DataFrame(stats).fillna(0).to_excel(out_writer, 'Stats')
    samples_out.to_excel(out_writer, casecon)
    matrix.to_excel(out_writer, 'Significant samples')
    matrix_pos.to_excel(out_writer, 'POS Significant samples')
    matrix_neg.to_excel(out_writer, 'NEG Significant samples')
    all_pwy_stats.to_excel(out_writer, 'ALL PWY STATS')


    out_writer.close()


if __name__=='__main__':
    main()
