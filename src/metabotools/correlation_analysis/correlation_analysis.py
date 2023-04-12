
import pandas as pd
from IPython import embed
from scipy import stats
import numpy as np
import os

def make_corr_spreadsheet(data, type, method, infile, study, db_file='', table1_file='', FDR_threshold=0.05):
    

    import datetime
    now = datetime.datetime.now()
    # embed()
    # if os.path.exists('tmp_corr.xlsx'):
    #     corrs = pd.read_excel('tmp_corr.xslx', index_col=0)
    # else:
    corrs = get_correlations(data, method=method)
    #     corrs.to_excel('tmp_corr.xlsx')
    

    DIR='\\'.join(__file__.split('\\')[:-1])+'\\'
    MTDB_DIR = '\\'.join(__file__.split('\\')[:-3])+'\\'
    if db_file=='':
        db_file = MTDB_DIR+'mtdb\\MasterDatabase_v5.1.2.xlsx'
    
    cyto_xy = DIR+'Cytoscape_xy_with_ID_v2.xlsx'
    
    mtdb_map = map_to_mdtb(corrs, db_file)

    cyto_map = map_to_cytoscape(mtdb_map, cyto_xy, method)['cyto']
    
    sig_corrs = corrs[corrs['FDR']<FDR_threshold]
    r_cutoff = sig_corrs['abs(r)'].min()
    neg_corr = sig_corrs[sig_corrs['r']<0]
    pos_corr = sig_corrs[sig_corrs['r']>0]
    
    version = '1.0'    

    readme = {'Input File':infile,
        'Date':str(now.isoformat()),
        'Type':type,
        'Correlation Method':method,
        'Metabolites Measured':len(data.columns),
        'Total Pairwise Correlations':len(corrs),
        'Sample Number':len(data),
        'Total Pairs <= FDR':len(sig_corrs),
        'r threshold at FDR = %0.2f'%FDR_threshold:'%0.2f'%r_cutoff,
        'Positively correlated pairs':len(pos_corr),
        'Negatively correlated pairs':len(neg_corr),
        'VRSC version':version,
        }
    # embed()
    
    if table1_file!='':
        readme['Z-Scores File']=table1_file
        table1 = pd.read_excel(table1_file,'Metabolites')
        alias = pd.read_excel(db_file,'Alias')
        mrm = pd.read_excel(db_file,'MRM')
        chemical = pd.read_excel(db_file,'Chemical')
        mtdb = pd.merge(alias, mrm)
        mtdb = pd.merge(mtdb,chemical)

        table1 = pd.merge(mtdb, table1)[['MRM Name','MRM Index','Z Score','Chemical Index','Cytoscape Name']]
        # embed()
        # new code:
        unity = pd.merge(table1[['Z Score','Cytoscape Name']], cyto_map, left_on='Cytoscape Name', right_on='name1', how='right')
        del unity['Cytoscape Name']
        unity = unity.rename(columns={'Z Score':'ZScore1'})
        unity2 = pd.merge(table1, unity, left_on='Cytoscape Name', right_on='name2', how='right')
        del unity2['Cytoscape Name']
        unity2 = unity2.rename(columns={'Z Score':'ZScore2'})
        

        unity_for_nate = unity2[['Pathway','name1','x1','y1','Chemical Index1','ZScore1','name2','x2','y2','Chemical Index2','ZScore2','Duplicate Pairs']]
        
        #old code:
        # unity = pd.merge(table1, cyto_map[['name1','x1','y1']], left_on='Cytoscape Name', right_on='name1').drop_duplicates()
        # unity = unity[['Chemical Index', 'Cytoscape Name', 'Z Score', 'x1','y1']]
        
    
    readme = pd.DataFrame(pd.Series(readme), columns=[''])
    
    writer = pd.ExcelWriter('%s_%s_%s_VR_SHEET_v%s.xlsx'%(study,type,method, version))
    
    readme.to_excel(writer, 'README')
    data.to_excel(writer,'Z-Scores')
    mtdb_map.to_excel(writer,'Pairwise Correlations')
    cyto_map.to_excel(writer, 'Cytoscape XY Map')
    if table1_file!='':
        unity2.to_excel(writer, 'Unity')
    writer.save()
    unity2.to_csv('%s_%s_%s_VR_CSV_FORJON.csv'%(study,type,method), sep=';')
    unity_for_nate.to_csv('%s_%s_%s_VR_CSV.csv'%(study,type,method))
    
    
def get_corrs(size, num_bootstraps, data, tables1, control, exp, type, db_file, method, Q_threshold, OUTDIR, header):
    '''
    This is called for bootstrapping analysis
    '''
    from itertools import combinations
    import random
    print('running %i bootstraps for size %i'%(num_bootstraps, size))
    out = pd.DataFrame()
    pos_out = pd.DataFrame()
    neg_out = pd.DataFrame()
    bootstrap_samples = pd.DataFrame()
    
    samples = data[data['Group']==exp].index.tolist()
    controls = data[data['Group']==control].index.tolist()
    
    stats = []
    
    controls = data[data['Group']==control].index.tolist()
    cases = data[data['Group']==exp].index.tolist()
    
    # determine number of possible combinations
    
    if type=='CONTROLS':
        print('Running for CONTROLS (n=%i)'%len(controls))
        if len(controls)-size<5:
            comb = combinations(controls, size)
            comb = [c for c in comb]
            comb.insert(0,controls) # add full set to first pos
        else:
            comb = [1]*1000
    elif type=='CASES':
        print('Running for CASES (n=%i)'%len(cases))
        if len(cases)-size<5: # number of combinations will be huge so just go forward
            comb = combinations(samples, size)
            comb = [c for c in comb]
            comb.insert(0,cases)
        else:
            comb = [1]*1000
    print('there are >%i unique combinations for this sample size'%len(comb))
    # embed()
    short=False
    if len(comb) < num_bootstraps:
        print('run fewer than max because # combos = %i < max boostrap (%i)'%(len(comb),num_bootstraps))
        print('thus, running %i bootstraps'%len(comb))
        num_bootstraps = len(comb)
        short = True # denote that we can run fewer than max because # combos < max boostrap
        
    def check_dupe_subset(subset, cur_samples):
        for c in cur_samples:
            if len(subset)==len(set(subset)&c):
                print('subset has already been run, drawing again randomly')
                return True # denotes that subset is a dupe
        else:
            return False
        
    
    all_pwy_stats = []
    cur_samples = []
    for i in range(0,num_bootstraps+1):
        print('size %i: bootstrap %i'%(size,i))
        
        # to check if samples have already been run:
        # TODO: ensure we don't get duplicate samples
        
        if i>1:
            for b in bootstrap_samples[0].unique():
                s = set(bootstrap_samples[bootstrap_samples[0]==b].index.tolist())
                cur_samples.append(s)
            
        if type=='CONTROLS':
            if short==True:
                subset = [c for c in comb[i-1]]
            else:
                dupe=True
                while dupe==True:
                    subset = random.sample(controls, size)
                    print('checking whether this subset has already been run')
                    dupe = check_dupe_subset(subset,cur_samples)
                # TODO ensure we dont get duplicate samples: not yet implemented - need to be done for cases below as well
                # for s in cur_samples:
                    # print(len(set(s)-set(subset)))
                    # subset = random.sample(controls, size) # resample if dupes
                
        elif type=='CASES':
            if short==True:
                subset = [c for c in comb[i-1]]
            else:
                dupe=True
                while dupe==True:
                    subset = random.sample(cases, size)
                    print('checking whether this subset has already been run')
                    dupe = check_dupe_subset(subset,cur_samples)
                
        # zscores for full dataaset
        # CALC Z SCORES BEFORE SUBSAMPLING:
        zsc = get_zscores(data, exp, control)
        zscores = pd.concat([zsc['case'],zsc['control']])
        zscores=pd.concat([data['Group'],zscores], axis=1)
        print('calculating z-scores once - PRE SUBSAMPLING')
        # embed()
        
        ## TMP DEBUG:: REMOVE ME! 3/17/2022
        # print('TEMP MARCH 16 - using raw AUCs rather than z-scores')
        # zscores = data
        
        
        if i==0:
            controls_subset = controls
        else:
            bootstrap_samples = pd.concat([bootstrap_samples, pd.DataFrame(['boostrap_%i'%i]*len(subset), subset)])
        
        # embed()
        
        if type=='CONTROLS':
            # data2 = data.loc[subset+samples,:] # subset of controls
            data2 = zscores.loc[subset+samples,:] # subset of controls
        elif type=='CASES':
            # data2 = data.loc[controls+subset,:] # subset of samples
            data2 = zscores.loc[controls+subset,:] # subset of samples
       
        print(len(data2))
        
        
        corrs = bootstrap(data2, control, exp, db_file, method=method, use_zscore=False) # this will now be TableS2
        
        if type=='CONTROLS':
            corr_results = corrs['%s_control'%method]
        elif type=='CASES':
            corr_results = corrs['%s_case'%method]
            
        # FILTERING:
        
        corr_results = corr_results[corr_results['Q value'] < Q_threshold]
        print('note: only filtering by Q values now TODO implement p and r val filtering')
       
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
            'std r':input['r'].std(),
            'mean r':input['r'].mean()
            }
            
            return stats
            
        cur_stats = get_stats(corr_results)
        cur_stats['type']='all'
        
        if i==0:
            cur_stats['bootstrap']='full_set'
        else:   
            cur_stats['bootstrap']='boostrap_%i'%i
        stats.append(cur_stats)
        
        # getting canonical pathway stats:
        for unique_pwy in tables1['Pathway Name'].unique():
            mets_in_pwy = tables1[tables1['Pathway Name']==unique_pwy]['MRM Name'].tolist()
            pwy_cors = corr_results[corr_results['from'].isin(mets_in_pwy)]            
            pwy_cors = pd.concat([pwy_cors,corr_results[corr_results['to'].isin(mets_in_pwy)]],axis=0)       
            pwy_stats = get_stats(pwy_cors)
            pwy_stats['pathway']=unique_pwy
            pwy_stats['bootstrap']='boostrap_%i'%i
            all_pwy_stats.append(pwy_stats)
        
        # getting knn pathway stats:
        for unique_pwy in tables1['kNN Cluster (of 10)'].unique():
            mets_in_pwy = tables1[tables1['kNN Cluster (of 10)']==unique_pwy]['MRM Name'].tolist()
            pwy_cors = corr_results[corr_results['from'].isin(mets_in_pwy)]            
            pwy_cors = pd.concat([pwy_cors, corr_results[corr_results['to'].isin(mets_in_pwy)]], axis=0)       
            pwy_stats = get_stats(pwy_cors)
            pwy_stats['pathway']='kNN cluster ' + str(unique_pwy)
            pwy_stats['bootstrap']='boostrap_%i'%i
            all_pwy_stats.append(pwy_stats)
            
        
        pos_corr_results = corr_results[corr_results['r'] >=0]
        neg_corr_results = corr_results[corr_results['r'] < 0]
        
        pos_stats = get_stats(pos_corr_results)
        pos_stats_in = get_stats(pos_corr_results[pos_corr_results['Pathway']=='IN'])
        pos_stats_out = get_stats(pos_corr_results[pos_corr_results['Pathway']=='OUT'])
        pos_stats['type']='positive'
        pos_stats_in['type']='positive'
        pos_stats_out['type']='positive'
        pos_stats['ptype']='ALL'
        pos_stats_in['ptype']='IN'
        pos_stats_out['ptype']='OUT'
        
        neg_stats = get_stats(neg_corr_results)
        neg_stats_in = get_stats(neg_corr_results[neg_corr_results['Pathway']=='IN'])
        neg_stats_out = get_stats(neg_corr_results[neg_corr_results['Pathway']=='OUT'])
        neg_stats['type']='negative'
        neg_stats_in['type']='negative'
        neg_stats_out['type']='negative'
        neg_stats['ptype']='ALL'
        neg_stats_in['ptype']='IN'
        neg_stats_out['ptype']='OUT'
        
        if i==0:
            pos_stats['bootstrap']='full_set'
            neg_stats['bootstrap']='full_set'
            pos_stats_in['bootstrap']='full_set'
            pos_stats_out['bootstrap']='full_set'
            neg_stats_in['bootstrap']='full_set'
            neg_stats_out['bootstrap']='full_set'
        else:   
            pos_stats['bootstrap']='boostrap_%i'%i
            neg_stats['bootstrap']='boostrap_%i'%i
            pos_stats_in['bootstrap']='boostrap_%i'%i
            pos_stats_out['bootstrap']='boostrap_%i'%i
            neg_stats_in['bootstrap']='boostrap_%i'%i
            neg_stats_out['bootstrap']='boostrap_%i'%i
        stats.append(pos_stats)
        stats.append(neg_stats)
        stats.append(pos_stats_in)
        stats.append(pos_stats_out)
        stats.append(neg_stats_in)
        stats.append(neg_stats_out)
        
        
        if i==0:
            corr_results = corr_results.assign(bootstrap='full_set')
            pos_corr_results = pos_corr_results.assign(bootstrap='full_set')
            neg_corr_results = neg_corr_results.assign(bootstrap='full_set')
            
        else:
            corr_results = corr_results.assign(bootstrap='boostrap_%i'%i)
            pos_corr_results = pos_corr_results.assign(bootstrap='boostrap_%i'%i)
            neg_corr_results = neg_corr_results.assign(bootstrap='boostrap_%i'%i)
            
            
        out = pd.concat([out,corr_results])
        pos_out = pd.concat([pos_out, pos_corr_results])
        neg_out = pd.concat([neg_out, neg_corr_results])
        
        print('\n\n ----------------------- \n')
        
    
        
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
    samples_out = pd.pivot_table(bootstrap_samples_out, index='index', columns='bootstrap', values='pres')

    all_pwy_stats = pd.DataFrame(all_pwy_stats)


    out_writer = pd.ExcelWriter(OUTDIR+'/%s_%i_bootstraps_%i_samples_%s.xlsx'%(type,num_bootstraps,size,method))
    pd.DataFrame(header).to_excel(out_writer, 'Header')
    pd.DataFrame(stats).to_excel(out_writer, 'Stats')
    samples_out.to_excel(out_writer, type)
    matrix.to_excel(out_writer, 'Significant samples')
    matrix_pos.to_excel(out_writer, 'POS Significant samples')
    matrix_neg.to_excel(out_writer, 'NEG Significant samples')
    all_pwy_stats.to_excel(out_writer, 'ALL PWY STATS')


    out_writer.close()
    
def get_upper_diag(df):
    # data = np.triu(np.ones(df.shape),k=1).astype(np.bool)
    df2 = df.where(np.triu(np.ones(df.shape),k=1).astype(np.bool))
    df3 = df2.stack().reset_index()
    df3.columns = ['Row','Column','Value']
    return df3

def run_tableS3(tableS2_file, case_con, filename='test.xlsx'):
    print(f'Running Table S3 on {tableS2_file} for {case_con}')
    data = pd.read_excel(tableS2_file,'Fishers Exact p - Pearson Out')
    data2 = pd.read_excel(tableS2_file,'pearson_'+case_con)
    ts2readme = pd.read_excel(tableS2_file,'README', index_col=0)
    
    data_out = data2[data2['Pathway']=='OUT']  # focus on out of pathway corrs
    stats05 = calc_ts3_stats(data, data_out, 0.05, case_con)
    stats0005 = calc_ts3_stats(data, data_out, 0.0005, case_con)
    
    qvq = pd.merge(stats05, stats0005)
    qvq['Node Robustness = (n at q<0.0005) / (n at q<0.05)'] = qvq[f'Nodes (metabolites) in the {case_con} risk hub (q<0.0005)'] / qvq[f'Nodes (metabolites) in the {case_con} risk hub (q<0.0500)']
    qvq['Edge Robustness = (n at q<0.0005) / (n at q<0.05)'] = qvq[f'Edges (correlations ) from the {case_con} risk hub (q<0.0005)'] / qvq[f'Edges (correlations ) from the {case_con} risk hub (q<0.0500)']

    qtec = get_q_thresh_edge_counts(data2)
    data.index = data['Pathway Name']
    qtec['Total Metabolites in Pathway'] = data['Number of chemicals measured in the pathway']

    writer = pd.ExcelWriter(filename)

    ts2readme.to_excel(writer, 'README')
    qvq.to_excel(writer, 'q 0.05 vs 0.0005')
    qtec.fillna(0).to_excel(writer, 'q threshold edge counts')

    writer.close()
    

def get_q_thresh_edge_counts(data2):
    '''
    Calculate significant edges at various q thresholds
    '''
    data_out = data2[data2['Pathway']=='OUT']  # focus on out of pathway corrs
    pos_out = data_out[data_out['r']>0]
    neg_out = data_out[data_out['r']<0]

    results = pd.DataFrame()
    thresholds = [0.05, 0.01, 0.005, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001]
    for threshold in thresholds:
        op = data_out[data_out['Q value']<threshold]
        opp = pos_out[pos_out['Q value']<threshold]
        opn = neg_out[neg_out['Q value']<threshold]

        pwy1 = op.groupby('Pathway1').count()['from']
        pwy2 = op.groupby('Pathway2').count()['from']
        op_total = pd.concat([pwy1, pwy2], axis=1).transpose().sum()
        op_total = pd.DataFrame(op_total)
        op_total['q'] = threshold
        op_total['subset'] = 'All Out'
        
        pwy1 = opp.groupby('Pathway1').count()['from']
        pwy2 = opp.groupby('Pathway2').count()['from']
        opp_total = pd.concat([pwy1, pwy2], axis=1).transpose().sum()
        opp_total = pd.DataFrame(opp_total)
        opp_total['q'] = threshold
        opp_total['subset'] = 'All Out +r'

        pwy1 = opn.groupby('Pathway1').count()['from']
        pwy2 = opn.groupby('Pathway2').count()['from']
        opn_total = pd.concat([pwy1, pwy2], axis=1).transpose().sum()
        opn_total = pd.DataFrame(opn_total)
        opn_total['q'] = threshold
        opn_total['subset'] = 'All Out -r'
        if len(opn_total)==0: # all are zero -> fill as zero
            opn_total = opp_total.copy()
            opn_total[0] = 0
            opn_total['subset'] = 'All Out -r'
        
        results = pd.concat([results, op_total, opp_total, opn_total])


    results = results.reset_index()

    resmat = pd.pivot_table(results, index='index', columns=['subset','q'], values=0).fillna(0)
    
    return resmat

def calc_ts3_stats(data, data2, q_value, case_con):


    data2 = data2[data2['Q value']<q_value]

    out = []

    for i in data.index:

        pwy = data.loc[i,'Pathway Name']
        meas_chem = data.loc[i,'Number of chemicals measured in the pathway']
        
        data3 = data2[data2['Pathway1']==pwy]
        data4 = data2[data2['Pathway2']==pwy]
        
        edges = len(data3)+len(data4)
        
        from_mets = data3.groupby('from').count().index.tolist()
        to_mets = data4.groupby('to').count().index.tolist()

        union = set(from_mets) | set(to_mets)
        
        out.append({
            'Pathway':pwy,
            'Total Metabolites in the pathway':meas_chem,
            f'Edges (correlations ) from the {case_con} risk hub (q<%0.4f)'%q_value:edges,
            f'Nodes (metabolites) in the {case_con} risk hub (q<%0.4f)'%q_value:len(union)
        })


    out = pd.DataFrame(out)    
    out['Connectivity = Edges/Nodes at q<%0.4f'%q_value] = out[f'Edges (correlations ) from the {case_con} risk hub (q<%0.4f)'%q_value]/out[f'Nodes (metabolites) in the {case_con} risk hub (q<%0.4f)'%q_value]
    
    # out.to_csv('table3_comparison.py')
        
    return out
    
def bootstrap(data, control, exp, db_file, log_transform=True, use_zscore=True, method='spearman'):
    '''
    Expects data in to have interger indexes and 2 columns:
    1 for sample name: "Sample"
    2 for group: "Group"
    '''
    
    # data = pd.read_csv(data_file)
    
    
    
    if use_zscore==True:
        zscores = get_zscores(data, exp, control, log_transform=log_transform)
        zscores = pd.concat([zscores['case'],zscores['control']], axis=0)
        zscores['Group']=data['Group']
        data = zscores
    else:
        data = data[data.columns[1:].tolist()+[data.columns[0]]] # if not calcing z-scores need to add 'Group' as last column
    # embed()
    
    targeted = len(data.columns) # minus 2 to account for Sample and Group columns
    
    counts = data.count()
    to_compare = counts[counts>0]
    data = data[to_compare.index]
    
    measured = len(data.columns) # minus 2 to account for Sample and Group columns
    
    # methods=['pearson','spearman']
    # methods=['spearman']
    
    
    case_data = data[data['Group']==exp]
    control_data = data[data['Group']==control]
    summaries = {}
    raw_data = {}
    
    expected = ((measured**2)-measured)/2
    
    # for method in methods:
    case_corr = get_correlations_w_missing_data(case_data, method=method)
    case_corr_w_pathways = map_to_mdtb(case_corr, db_file)
    if len(case_corr)!=len(case_corr_w_pathways):
        print('error merging case corr list with MDTB')
        embed()
    
    raw_data[method+'_case'] = case_corr_w_pathways
    # case_results = summarize_results(case_corr_w_pathways, 'Case', method)
    # case_totals = [c for c in case_results.columns if 'Total' in c]
    # case_total_sum = case_results[case_totals].transpose().sum()
    
    control_corr = get_correlations_w_missing_data(control_data, method=method)
    control_corr_w_pathways = map_to_mdtb(control_corr, db_file)
    if len(control_corr)!=len(control_corr_w_pathways):
        print('error merging contol corr list with MDTB')
        embed()
    raw_data[method+'_control'] = control_corr_w_pathways
    
    # control_results = summarize_results(control_corr_w_pathways, 'Control', method)
    # control_totals = [c for c in control_results.columns if 'Total' in c]
    # control_total_sum = control_results[control_totals].transpose().sum()
    
    # summary = pd.merge(case_results, control_results, left_index=True, right_index=True)
    # diff = case_total_sum-control_total_sum
    # summary['Difference in Total Correlations (q<0.05) by %s'%method] = diff
    
    
    # summaries[method] = summary
    
    return raw_data
        
        # embed()
        
    # corr_summary = pd.merge(summaries[methods[0]], summaries[methods[1]], left_index=True, right_index=True)
    # corr_summary = corr_summary.fillna(0)
    
    # return corr_summary
    
    # embed()
    
    
def correlation_summary_tables2(data_file, control, exp, db_file, outfile_name, log_transform=True, use_zscore=True, q_filter=0.05, r_filter=0):
    '''
    Expects data in to have interger indexes and 2 columns:
    1 for sample name: "Sample"
    2 for group: "Group"
    '''
    
    # data = pd.read_csv(data_file)
    
    data = pd.read_csv(data_file, index_col=0)
    
    if use_zscore==True:
        zscores = get_zscores(data, exp, control, log_transform=log_transform)
        zscores = pd.concat([zscores['case'], zscores['control']])
        zscores['Group']=data['Group']
        data = zscores
    else:
        data = data[data.columns[1:].tolist()+[data.columns[0]]] # if not calcing z-scores need to add 'Group' as last column
    # embed()
    
    targeted = len(data.columns)-1 # minus 2 to account for Sample and Group columns
    
    counts = data.count()
    to_compare = counts[counts>0]
    data = data[to_compare.index]
    
    measured = len(data.columns)-1 # minus 2 to account for Sample and Group columns

    
    methods=['pearson','spearman']
    
    
    case_data = data[data['Group']==exp]
    control_data = data[data['Group']==control]
    summaries = {}
    raw_data = {}
    
    expected = ((measured**2)-measured)/2
    
    for method in methods:
        case_corr = get_correlations_w_missing_data(case_data, method=method)
        case_corr_w_pathways = map_to_mdtb(case_corr, db_file)

        if len(case_corr)!=len(case_corr_w_pathways):
            print('error merging case corr list with MDTB')
            embed()
        
        raw_data[method+'_case'] = case_corr_w_pathways
        case_results = summarize_results(case_corr_w_pathways, 'Case', method, q_filter=q_filter, r_filter=r_filter)
        case_totals = [c for c in case_results.columns if 'r in Case Network' in c and 'Number of' in c]
        case_total_sum = case_results[case_totals].transpose().sum()
        
        control_corr = get_correlations_w_missing_data(control_data, method=method)
        control_corr_w_pathways = map_to_mdtb(control_corr, db_file)
        if len(control_corr)!=len(control_corr_w_pathways):
            print('error merging contol corr list with MDTB')
            embed()
        raw_data[method+'_control'] = control_corr_w_pathways
        
        control_results = summarize_results(control_corr_w_pathways, 'Control', method, q_filter=q_filter, r_filter=r_filter)
        control_totals = [c for c in control_results.columns if 'r in Control Network' in c and 'Number of' in c]
        control_total_sum = control_results[control_totals].transpose().sum()

        summary = pd.merge(case_results, control_results, left_index=True, right_index=True)
        diff = case_total_sum-control_total_sum
        summary['Difference in Total Correlations (q<%0.4f) by %s'%(q_filter, method)] = diff

        control_out_totals = [c for c in control_results.columns if 'Number of "Out of pathway" correlations in' in c][0]
        case_out_totals = [c for c in case_results.columns if 'Number of "Out of pathway" correlations' in c][0]
        out_diff = case_results[case_out_totals] - control_results[control_out_totals]
        summary['Difference in Out of Pathway Correlations (q<%0.4f) by %s'%(q_filter, method)] = out_diff

        control_in_totals = [c for c in control_results.columns if 'Number of "In-pathway" correlations in' in c][0]
        case_in_totals = [c for c in case_results.columns if 'Number of "In-pathway" correlations' in c][0]
        in_diff = case_results[case_in_totals] - control_results[control_in_totals]
        summary['Difference in In Pathway Correlations (q<%0.4f) by %s'%(q_filter, method)] = in_diff

        
        
        
        summaries[method] = summary

    # embed()
        

    corr_summary = pd.merge(summaries[methods[0]], summaries[methods[1]], left_index=True, right_index=True)
    corr_summary = corr_summary.fillna(0)

    alias = pd.read_excel(db_file,'Alias')
    mrm = pd.read_excel(db_file,'MRM')
    chemical = pd.read_excel(db_file,'Chemical', index_col=0)
    pathway = pd.read_excel(db_file,'Pathway')

    mtdb = pd.merge(alias, mrm, left_on='MRM Index', right_on='MRM Index')
    # embed()
    mtdb = pd.merge(mtdb, chemical, left_on='Chemical Index', right_on='Chemical Index')
    mtdb = pd.merge(mtdb, pathway, left_on='Pathway Index', right_on='Pathway Index')

    corr_summary = pd.merge(corr_summary, mtdb[['MRM Name_x','Pathway Name']], left_index=True, right_on='MRM Name_x')
    corr_summary.index = corr_summary['MRM Name_x']
    corr_summary.index.name = 'MRM Name'
    corr_summary.drop(columns=['MRM Name_x'], inplace=True)
    met_counts = corr_summary.groupby(['Pathway Name']).count()[corr_summary.columns[0]]

    for corr_type in ['spearman','pearson']:
        column_name = 'Difference in Total Correlations (q<%0.4f) by %s'%(q_filter, corr_type)
        dif_sums = corr_summary[[column_name,'Pathway Name']]
        dfs_pos = dif_sums[dif_sums[column_name]>0]
        dfs_neg = dif_sums[dif_sums[column_name]<0]
        
        dfs_pos = dfs_pos.groupby('Pathway Name').sum()
        norm_dfs_pos = dfs_pos['Difference in Total Correlations (q<%0.4f) by %s'%(q_filter,corr_type)]/met_counts
        dfs_pos = pd.concat([dfs_pos, norm_dfs_pos], axis=1)
        dfs_pos = dfs_pos.rename(columns={column_name:'Sum of positive differences - %s'%corr_type, 0:'NORMALIZED Sum of positive differences - %s'%corr_type})
        dfs_pos = dfs_pos.fillna(0)
        corr_summary = pd.merge(corr_summary, dfs_pos, left_on='Pathway Name', right_index=True, how='left')
        
        dfs_neg = dfs_neg.groupby('Pathway Name').sum().abs()
        norm_dfs_neg = dfs_neg['Difference in Total Correlations (q<%0.4f) by %s'%(q_filter,corr_type)]/met_counts
        dfs_neg = pd.concat([dfs_neg, norm_dfs_neg], axis=1)
        dfs_neg = dfs_neg.rename(columns={column_name:'Sum of abs(negative) differences - %s'%corr_type, 0:'NORMALIZED Sum of negative differences - %s'%corr_type})
        dfs_neg = dfs_neg.fillna(0)
        corr_summary = pd.merge(corr_summary, dfs_neg, left_on='Pathway Name', right_index=True, how='left')

    # embed()
    # add fishers exact calcs - 09212022 only for pearson (for now)
    pos_case = [c for c in corr_summary.columns if c.startswith('Number of (+) r for "Out of pathway" in Case Network pearson Correlations')][0]
    neg_case = [c for c in corr_summary.columns if c.startswith('Number of (-) r for "Out of pathway" in Case Network pearson Correlations')][0]
    pos_control = [c for c in corr_summary.columns if c.startswith('Number of (+) r for "Out of pathway" in Control Network pearson Correlations')][0]
    neg_control = [c for c in corr_summary.columns if c.startswith('Number of (-) r for "Out of pathway" in Control Network pearson Correlations')][0]
    case_all = [c for c in corr_summary.columns if c.startswith('Number of "Out of pathway" correlations in Case Network by pearson')][0]
    control_all = [c for c in corr_summary.columns if c.startswith('Number of "Out of pathway" correlations in Control Network by pearson')][0]
    out_pwy_diff = [c for c in corr_summary.columns if c.startswith('Difference in Out of Pathway Correlations')][0]

    cols = [case_all, control_all, pos_case, neg_case, pos_control, neg_control, out_pwy_diff]
    

    
    diff_col = 'Difference in Total Correlations (q<%0.4f) by pearson'%q_filter
    for_fisher = corr_summary.drop_duplicates().groupby('Pathway Name').sum()[cols]
    
    from scipy.stats import fisher_exact
    fisher_out = []
    for i in for_fisher.index:
        tmp = for_fisher.loc[i,:]
        table = [[tmp.loc[pos_case],tmp.loc[neg_case]], [tmp.loc[pos_control],tmp.loc[neg_control]]]
        oddsr, p = fisher_exact(table)
        
        corr_sum_pwy = corr_summary[corr_summary['Pathway Name']==i]
        mgi = corr_sum_pwy[diff_col].max()
        mgi_met = corr_sum_pwy[diff_col].idxmax()
        mdi = corr_sum_pwy[diff_col].min()
        mdi_met = corr_sum_pwy[diff_col].idxmin()
        fisher_out.append({
            'Pathway Name':i,
            'Network Quality by +/- p (Fishers Exact)':p,
            'Metabolite in the Pathway with the greatest INCREASE in out of pathway correlations difference (+ and - r; ASD-TD)  (q<%0.4f Pearson)'%q_filter:mgi_met,
            'Number of ASD - TD  Pearson Differences for increased metabolite':mgi,
            'Metabolite in the Pathway with the greatest DECREASE in out of pathway correlations difference (+ and - r; ASD-TD)  (q<%0.4f Pearson)'%q_filter:mdi_met,
            'Number of ASD - TD  Pearson Differences for decreased metabolite':mdi

        })

        
    fisher_out = pd.DataFrame(fisher_out)
    fisher_out = pd.merge(for_fisher, fisher_out, left_index=True, right_on='Pathway Name')

    mets_per_pwy = corr_summary.groupby('Pathway Name').count()[corr_summary.columns[0:2]]
    mets_per_pwy = mets_per_pwy.rename(columns={corr_summary.columns[0]:'Number of chemicals measured in the pathway'})
    fisher_out = pd.merge(mets_per_pwy, fisher_out, left_index=True, right_on='Pathway Name')
    fisher_out['Number of Out of pathway metabolites'] = measured-fisher_out['Number of chemicals measured in the pathway']
    fisher_out['Total number of possible out of pathway correlations'] = fisher_out['Number of chemicals measured in the pathway']*fisher_out['Number of Out of pathway metabolites']

    fisher_out['Out of Pathway Difference (Case-Control)'] = fisher_out[out_pwy_diff]
    fisher_out['Difference in positive (+r) correl (Case-Control)'] = fisher_out[pos_case]-fisher_out[pos_control]
    fisher_out['Difference in negative (-r) correl (Case-Control)'] = fisher_out[neg_case]-fisher_out[neg_control]

    cols_to_keep = [
        'Pathway Name', 
        'Network Quality by +/- p (Fishers Exact)',
        'Number of chemicals measured in the pathway',
        'Number of Out of pathway metabolites',
        'Total number of possible out of pathway correlations',
        case_all,
        control_all,
        pos_case,
        neg_case,
        pos_control,
        neg_control,
        out_pwy_diff,
        'Difference in positive (+r) correl (Case-Control)',
        'Difference in negative (-r) correl (Case-Control)',
        'Metabolite in the Pathway with the greatest INCREASE in out of pathway correlations difference (+ and - r; ASD-TD)  (q<%0.4f Pearson)'%q_filter,
        'Number of ASD - TD  Pearson Differences for increased metabolite',
        'Metabolite in the Pathway with the greatest DECREASE in out of pathway correlations difference (+ and - r; ASD-TD)  (q<%0.4f Pearson)'%q_filter,
        'Number of ASD - TD  Pearson Differences for decreased metabolite',
        
        
    ]

    fisher_out = fisher_out[cols_to_keep]

    # ensure that case_all = case_pos + case_neg, same for control
    check = (corr_summary[pos_case] + corr_summary[neg_case]) - corr_summary[case_all]
    if check.sum()!=0:
        print('error, case pos and neg != all case')
        embed()
    check = (corr_summary[pos_control] + corr_summary[neg_control]) - corr_summary[control_all]
    if check.sum()!=0:
        print('error, case pos and neg != all control')
        embed()
    

    venn_res = venn_analysis(raw_data['pearson_case'], raw_data['pearson_control'], q_filter, r_filter)
    

    readme = {'README':{
        'Unique correlations':len(case_corr),
        'expected':expected,
        'infile':data_file,
        'measured':measured,
        'targeted':targeted,
        'Data format analyzed':'Zlogs',
        'N (case=%s)'%exp:len(case_data),
        'N (control=%s)'%control:len(control_data),
        'Notes':'',
        'Summary q filter <':q_filter,
        'Summary r filter abs(r)<':r_filter,
        'OOP Edges':len(raw_data['pearson_case'][raw_data['pearson_case']['Pathway']=='OUT'])
    }}

    readme = pd.DataFrame(readme)

    

    # naviaux additions 1/30/2023

    # changes to Correlations Summary tab:
    corr_summary2 = corr_summary.reset_index().drop_duplicates()
    corr_summary2 = pd.merge(corr_summary2, fisher_out[['Pathway Name','Number of chemicals measured in the pathway']])
    corr_summary2['Number of Out of pathway metabolites'] = len(corr_summary2) - corr_summary2['Number of chemicals measured in the pathway']
    corr_summary2['Case Number of non-significant Out Edges by %s (q>%0.4f)'%('pearson', q_filter)] = corr_summary2['Number of Out of pathway metabolites'] - corr_summary2['Number of "Out of pathway" correlations in %s Network by %s (q<%0.4f)'%('Case', 'pearson', q_filter)]
    corr_summary2['Control Number of non-significant Out Edges by %s (q>%0.4f)'%('pearson', q_filter)] = corr_summary2['Number of Out of pathway metabolites'] - corr_summary2['Number of "Out of pathway" correlations in %s Network by %s (q<%0.4f)'%('Control', 'pearson', q_filter)]
    # embed()
    from scipy.stats import fisher_exact
    
    total_corr = expected
    # embed()
    for i in corr_summary2.index:
            
        #Network Quality-- (Fishers Exact p) [Use columns K&L vs V&W]]
        table = [
                    [
                        corr_summary2.loc[i,'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)],
                        corr_summary2.loc[i,'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)]
                    ],
                    [
                        corr_summary2.loc[i,'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)],
                        corr_summary2.loc[i,'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quality = p
        
        # Network Quantity-- (Fishers Exact p) [Use columns I&J vs T&U]
        table = [
                    [
                        corr_summary2.loc[i,'Number of "Out of pathway" correlations in %s Network by %s (q<%0.4f)'%('Case', 'pearson', q_filter)],
                        corr_summary2.loc[i,'Case Number of non-significant Out Edges by %s (q>%0.4f)'%('pearson', q_filter)]
                    ],
                    [
                        corr_summary2.loc[i,'Number of "Out of pathway" correlations in %s Network by %s (q<%0.4f)'%('Control', 'pearson', q_filter)],
                        corr_summary2.loc[i,'Control Number of non-significant Out Edges by %s (q>%0.4f)'%('pearson', q_filter)]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quantity = p
        
        corr_summary2.loc[i,'Network Quality-- (Fishers Exact p)'] = network_quality
        corr_summary2.loc[i,'Network Quantity-- (Fishers Exact p)'] = network_quantity

    # jan 27 changes:
    corr_summary2['Number of Non-significant (+) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f)'%(q_filter)] = corr_summary2['Number of Out of pathway metabolites'] - corr_summary2['Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)]
    corr_summary2['Number of Non-significant (-) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f)'%(q_filter)] = corr_summary2['Number of Out of pathway metabolites'] - corr_summary2['Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)]
    corr_summary2['Number of Non-significant (+) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f)'%(q_filter)] = corr_summary2['Number of Out of pathway metabolites'] - corr_summary2['Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)]
    corr_summary2['Number of Non-significant (-) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f)'%(q_filter)] = corr_summary2['Number of Out of pathway metabolites'] - corr_summary2['Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)]


    for i in corr_summary2.index:
            
        #Network Quality-- (Fishers Exact p) 
        table = [
                    [
                        corr_summary2.loc[i,'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)],
                        corr_summary2.loc[i,'Number of Non-significant (+) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f)'%q_filter]
                    ],
                    [
                        corr_summary2.loc[i,'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)],
                        corr_summary2.loc[i,'Number of Non-significant (+) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f)'%q_filter]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quality_pos = p
        
        # Network Quantity-- (Fishers Exact p) 
        table = [
                    [
                        corr_summary2.loc[i,'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)],
                        corr_summary2.loc[i,'Number of Non-significant (-) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f)'%q_filter]
                    ],
                    [
                        corr_summary2.loc[i,'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)],
                        corr_summary2.loc[i,'Number of Non-significant (-) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f)'%q_filter]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quality_neg = p
        
        corr_summary2.loc[i,'Network Quality-- +r only, Sig v Non-sig (Fishers Exact p)'] = network_quality_pos
        corr_summary2.loc[i,'Network Quality-- (-)r only, Sig v Non-sig (Fishers Exact p)'] = network_quality_neg
    

    # changes to Fishers Exact p - Pearson Out tab:
    fisher_out['Number of Non-Significant "Out of pathway" correlations in Case Network by %s (q>%0.4f)'%('pearson', q_filter)] = fisher_out['Total number of possible out of pathway correlations'] - fisher_out['Number of "Out of pathway" correlations in Case Network by %s (q<%0.4f)'%('pearson', q_filter)]
    fisher_out['Number of Non-Significant "Out of pathway" correlations in Control Network by %s (q>%0.4f)'%('pearson', q_filter)] = fisher_out['Total number of possible out of pathway correlations'] - fisher_out['Number of "Out of pathway" correlations in Control Network by %s (q<%0.4f)'%('pearson', q_filter)]
    for i in fisher_out.index:
            
        #Network Quality-- (Fishers Exact p) [Use columns K&L vs V&W]]
        table = [
                    [
                        fisher_out.loc[i,'Number of "Out of pathway" correlations in Case Network by %s (q<%0.4f)'%('pearson', q_filter)],
                        fisher_out.loc[i,'Number of Non-Significant "Out of pathway" correlations in Case Network by %s (q>%0.4f)'%('pearson', q_filter)]
                    ],
                    [
                        fisher_out.loc[i,'Number of "Out of pathway" correlations in Control Network by %s (q<%0.4f)'%('pearson', q_filter)],
                        fisher_out.loc[i,'Number of Non-Significant "Out of pathway" correlations in Control Network by %s (q>%0.4f)'%('pearson', q_filter)]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quantity = p
        print(p)
        fisher_out.loc[i,'Network Quantity-- (Fishers Exact p)'] = network_quantity

    # jan 27 changes:


    fisher_out['Number of Non-significant (+) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)] = fisher_out['Total number of possible out of pathway correlations'] - fisher_out['Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)]
    fisher_out['Number of Non-significant (-) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)] = fisher_out['Total number of possible out of pathway correlations'] - fisher_out['Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)]
    fisher_out['Number of Non-significant (+) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)] = fisher_out['Total number of possible out of pathway correlations'] - fisher_out['Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)]
    fisher_out['Number of Non-significant (-) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)] = fisher_out['Total number of possible out of pathway correlations'] - fisher_out['Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)]

    for i in fisher_out.index:
            
        #Network Quality-- (Fishers Exact p) [Use columns K&L vs V&W]]
        table = [
                    [
                        fisher_out.loc[i,'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)],
                        fisher_out.loc[i,'Number of Non-significant (+) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)]
                    ],
                    [
                        fisher_out.loc[i,'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)],
                        fisher_out.loc[i,'Number of Non-significant (+) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quality = p
        
        # Network Quantity-- (Fishers Exact p) [Use columns I&J vs T&U]
        table = [
                    [
                        fisher_out.loc[i,'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Case', 'pearson', q_filter, total_corr)],
                        fisher_out.loc[i,'Number of Non-significant (-) r for "Out of pathway" in Case Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)]
                    ],
                    [
                        fisher_out.loc[i,'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%('Control', 'pearson', q_filter, total_corr)],
                        fisher_out.loc[i,'Number of Non-significant (-) r for "Out of pathway" in Control Network pearson Correlations (q>%0.4f) of %i Total'%(q_filter, total_corr)]
                    ]
                ] 
        oddsr, p = fisher_exact(table)
        network_quantity = p
        
        fisher_out.loc[i,'Network Quality-- +r only, Sig v Non-sig (Fishers Exact p)'] = network_quality
        fisher_out.loc[i,'Network Quality-- (-)r only, Sig v Non-sig (Fishers Exact p)'] = network_quantity


    writer = pd.ExcelWriter(outfile_name)
    readme.to_excel(writer,'README')
    for k in raw_data.keys():
        if expected==len(raw_data[k]):
            raw_data[k].to_excel(writer,k)
        else:
            print('raw data for %s doesnt match expected number of corrs'%k)
            embed()

    corr_summary2.to_excel(writer,'Correlations Summary', index=False)

    ''' adding column Bob requested 2/9/23'''

    data = fisher_out
    data2 = raw_data['pearson_case']
    
    data2 = data2[data2['Q value']<0.05] # hardcoded here and below in label
    data2 = data2[data2['Pathway']=='OUT']

    for i in data.index:

        pwy = data.loc[i,'Pathway Name']
        
        data3 = data2[data2['Pathway1']==pwy]
        data4 = data2[data2['Pathway2']==pwy]
        
        from_mets = data3.groupby('from').count().index.tolist()
        to_mets = data4.groupby('to').count().index.tolist()

        union = set(from_mets) | set(to_mets)
        
        data.loc[i,'Nodes (metabolites) in the ASD risk hub (q<0.05)'] = len(union)
        
    # embed()

    ''' END BOB 2/9/23 addition '''  


    fisher_out.to_excel(writer,'Fishers Exact p - Pearson Out')
    if use_zscore==True:
        zscores.to_excel(writer,'Zscores')
    venn_res.to_excel(writer, 'venn_comparison_out_of_pwy')
    writer.close()
        
    print('saved results to %s'%outfile_name)
    
def venn_analysis(venn_case, venn_control, q_filter, r_filter):

    venn_case = venn_case[venn_case['Pathway']=='OUT']
    venn_control = venn_control[venn_control['Pathway']=='OUT']

    '''
    which are DIFFERENT between cases and controls. The list should have 5 categories:
a. Different in kind: edges present in cases but NOT in controls
b. Different in kind: edges absent in cases but present in controls
c. Different in sign: + in cases, but "-" in controls
d. Different in sign: - in cases, but "+" in controls
e. SAME in kind and sign
'''
    

    venn_case_sig = venn_case[venn_case['Q value']<q_filter]
    venn_case_sig['name'] = venn_case_sig['from'] + '<->' + venn_case_sig['to']
    venn_case_pos = venn_case_sig[venn_case_sig['r']>=r_filter]
    venn_case_neg = venn_case_sig[venn_case_sig['r']<=-1*r_filter]

    venn_control_sig = venn_control[venn_control['Q value']<q_filter]
    venn_control_sig['name'] = venn_control_sig['from'] + '<->' + venn_control_sig['to']
    venn_control_pos = venn_control_sig[venn_control_sig['r']>=r_filter]
    venn_control_neg = venn_control_sig[venn_control_sig['r']<=-1*r_filter]

    case = set(venn_case_sig['name'].tolist())
    control = set(venn_control_sig['name'].tolist())

    pos_case = set(venn_case_pos['name'].tolist())
    pos_control = set(venn_control_pos['name'].tolist())

    neg_case = set(venn_case_neg['name'].tolist())
    neg_control = set(venn_control_neg['name'].tolist())

    a = case-control
    b = control-case
    c = pos_case & neg_control
    d = neg_case & pos_control
    e = pos_case & pos_control
    f = neg_case & neg_control

    out = {
        'a. Different in kind: edges present in cases but NOT in controls':pd.Series(list(a)),
        'b. Different in kind: edges absent in cases but present in controls':pd.Series(list(b)),
        'c. Different in sign: + in cases, but "-" in controls':pd.Series(list(c)), 
        'd. Different in sign: - in cases, but "+" in controls':pd.Series(list(d)), 
        'e. SAME in kind and sign (pos)':pd.Series(list(e)), 
        'f. SAME in kind and sign (neg)':pd.Series(list(f))}

    out = pd.DataFrame(out)
    return out

    
def summarize_results(corr, sample, method, q_filter = 0.05, r_filter = 0):
    total_corr = len(corr)
    sig_by_q = corr[corr['Q value']<=q_filter] # q filter
    sig_by_q = sig_by_q[sig_by_q['abs(r)']>=r_filter] # r filter
    # neg_corr = sig_by_q[sig_by_q['r']<=-1*r_filter]

    unique_mets = set(corr['from'].unique().tolist() + corr['to'].unique().tolist())


    def combine_from_to(data, c):
        tmp1 = data[data['from']==c]
        tmp2 = data[data['to']==c]
        return pd.concat([tmp1,tmp2])
    out = {}
    all_mets = corr['from'].unique().tolist() + corr['to'].unique().tolist()
    all_mets = list(set(all_mets))
    for c in all_mets:
        sig_q = combine_from_to(sig_by_q,c)
        
        in_pwy = sig_q[sig_q['Pathway']=='IN']
        out_pwy = sig_q[sig_q['Pathway']=='OUT']

        pos_corr = sig_q[sig_q['r']>=r_filter]
        neg_corr = sig_q[sig_q['r']<=-1*r_filter]

        pos_corr_in = in_pwy[in_pwy['r']>=r_filter]
        pos_corr_out = out_pwy[out_pwy['r']>=r_filter]
        neg_corr_in = in_pwy[in_pwy['r']<=-1*r_filter]
        neg_corr_out = out_pwy[out_pwy['r']<=-1*r_filter]

        # pos = combine_from_to(pos_corr,c)
        # neg = combine_from_to(neg_corr,c)
        
        # pos_corr_in = pos_corr[pos_corr['Pathway']=='IN']
        # pos_corr_out = pos_corr[pos_corr['Pathway']=='OUT']
        # neg_corr_in = neg_corr[neg_corr['Pathway']=='IN']
        # neg_corr_out = neg_corr[neg_corr['Pathway']=='OUT']

        # pos_corr_in_cnt = combine_from_to(pos_corr_in,c)
        # pos_corr_out_cnt = combine_from_to(pos_corr_out,c)
        # neg_corr_in_cnt = combine_from_to(neg_corr_in,c)
        # neg_corr_out_cnt = combine_from_to(neg_corr_out,c)
        
        out[c] = {
            'Median of (+) r in %s Network by %s (q<%.4f) '%(sample, method, q_filter): pos_corr['r'].median(),
            'Median of (-) r in %s Network by %s (q<%0.4f) '%(sample, method, q_filter): neg_corr['r'].median(),
            'Number of (+) r in %s Network %s Correlations (q<%0.4f) of %i Total'%(sample, method, q_filter, total_corr): len(pos_corr),
            'Number of (-) r in %s Network %s Correlations (q<%0.4f) of %i Total'%(sample, method, q_filter, total_corr): len(neg_corr),
            'Number of "In-pathway" correlations in %s Network by %s (q<%0.4f) '%(sample, method, q_filter): len(in_pwy),
            'Number of "Out of pathway" correlations in %s Network by %s (q<%0.4f)'%(sample, method, q_filter): len(out_pwy),
            'Number of (+) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%(sample, method, q_filter, total_corr): len(pos_corr_out),
            'Number of (-) r for "Out of pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%(sample, method, q_filter, total_corr): len(neg_corr_out),
            'Number of (+) r for "In-pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%(sample, method, q_filter, total_corr): len(pos_corr_in),
            'Number of (-) r for "In-pathway" in %s Network %s Correlations (q<%0.4f) of %i Total'%(sample, method, q_filter, total_corr): len(neg_corr_in)
        }
    out = pd.DataFrame(out).transpose()
    if len(out)!= len(unique_mets):
        summarized_mets = out.index.tolist()
        missing_mets = unique_mets-set(summarized_mets)
        print(missing_mets, 'had no sig correlations')
        for m in list(missing_mets):
            out.loc[m,:]=np.nan
    
    return out
    
def get_correlations_w_missing_data(data, method='pearson'): # missing data slows this down significantly (filtering the NAs is slow)
    missing_data=True
    
    
    # BELOW WAS WRITTEN TO SPEED CORR CODE:
    corr = data.corr(method=method)
    
    if method=='pearson':
        pvals = data.corr(method=lambda x, y: stats.pearsonr(x, y)[1]) - np.eye(*corr.shape)
    elif method=='spearman':
        pvals = data.corr(method=lambda x, y: stats.spearmanr(x, y)[1]) - np.eye(*corr.shape)
    else:
        print('ERROR in get_correlations_w_missing_data(): invalid method')
        return
    
    
    corr_up = corr.where(np.triu(np.ones(corr.shape)).astype(np.bool))
    corr_up_stacked = corr_up.stack().reset_index()
    corr_up_stacked.columns = ['from','to','r']
    corr_up_stacked = corr_up_stacked[corr_up_stacked['from']!=corr_up_stacked['to']]

    pvals_up = pvals.where(np.triu(np.ones(pvals.shape)).astype(np.bool))
    pvals_up_stacked = pvals_up.stack().reset_index()
    pvals_up_stacked.columns = ['from','to','p']
    pvals_up_stacked = pvals_up_stacked[pvals_up_stacked['from']!=pvals_up_stacked['to']]

    out = pd.merge(corr_up_stacked, pvals_up_stacked)
    
    # OLD CORR CODE:
    '''
    upper = get_upper_diag(corr)
    
    out = []
    # idx=0
    # for i in corr.index:
        # for c in corr.columns:
    debug = pd.ExcelWriter('debug.xlsx')
    count=1
    for i in upper.index:
        row_i = upper.loc[i,'Row']
        col_i = upper.loc[i,'Column']
        # embed()
        # name = []
        # name.append(row_i)
        # name.append(col_i)
        # name.sort()
        # name_out='|'.join(name)
        
        if missing_data==True:
            
            row = data[row_i].dropna()
            col = data[col_i].dropna()
            
            shared = list(set(row.index)&set(col.index))
            # embed()
            row = row.loc[shared]
            col = col.loc[shared]
            # if c.startswith('T___'):
                # embed()
        else:
            row=data[row_i]
            col=data[col_i]
        
        if len(row)>3 and len(col)>3:
            if method=='pearson':
                r, p = stats.pearsonr(row, col)
            elif method=='spearman':
                r, p = stats.spearmanr(row, col)
            if r==-1:
                # embed()
                tmp = pd.DataFrame([row,col]).transpose()
                tmp.to_excel(debug,'%i'%(count))
                count+=1
                
            out.append({'from':row_i,'to':col_i, 'r':r, 'p':p})
            # out.append({'name':name_out, 'r':r, 'p':p})
        # idx+=1
    out = pd.DataFrame(out)
    
    '''
    # END OLD CORR CODE
    
    # debug.close()
    measured = len(data.columns)-1
    expect = (measured**2-measured)/2
    
    print( 'input cols', len(data.columns))
    print('list len:', len(out), 'expect', expect)
    '''
    out = out.drop_duplicates()
    print 'drop dupes:', len(out), 'expect', expect/2
    
    # NOTE drop dupes doesnt quite work for spearman - numerical imprecision seems to cause some dupes to remain:
    tmp = out['name'].drop_duplicates()
    out = out.ix[tmp.index]
    print 'drop dupes by name:', len(out)
    
    out['from']=out['name'].str.split('|',expand=True)[0]
    out['to']=out['name'].str.split('|',expand=True)[1]
    out = out[out['r']<=0.99999]
    print 'drop identity:', len(out), 'expect', expect/2-len(data.columns)
    
    del out['name']
    '''

    # FDR:
    from statsmodels.stats.multitest import fdrcorrection_twostage
    res = fdrcorrection_twostage(out['p'], alpha=0.05, method='bh')
    rej = res[0]
    pval_corr = res[1]

    # Qvalue:
    from .. import qvalue
    
    ix = out[out['p']>0.999999].index # check if any p-values=1 -> causes an error in qvalues
    if len(ix)>0:
        out.loc[ix,'p'] = 0.999999
    
    qv = qvalue.estimate(out['p']) # Storey's q value
    

    out['FDR'] = pval_corr
    out['Q value'] = qv
    out['abs(r)'] = abs(out['r'])
    # out.to_excel(debug)
    # debug.close()
    

    return out

def get_zscores_old(data, case, control):
    # Log transform data
    log2 = np.log2(data[data.columns[1:]])
    log2['Group']=data['Group']
    data = log2

    # split case and controls
    control = data[data['Group']==control][data.columns[:-1]]
    case = data[data['Group']==case][data.columns[:-1]]

    # calculate z-scores
    control_z = pd.DataFrame(stats.zscore(control, ddof=1), index=control.index, columns=control.columns)
    control_mean = control.mean()
    control_std = control.std()
    case_z = (case-control_mean)/control_std
    # embed()
    # only return cols with data:
    show=[]
    for c in data.columns[:-1]:
        if str(data[c].mean())!='nan':
            show.append(c)
    # embed()
            
    # control_z[show].to_csv('GWI_Control_zscore.csv')
    # case_z[show].to_csv('GWI_GWS_zscore.csv')

    
    return {'case':case_z[show], 'control':control_z[show]}

def get_zscores(data, case, control, log_transform=True):
        
    data2=pd.DataFrame()
    for c in data.columns[1:]:
        data2 = pd.concat([data2, pd.to_numeric(data[c])], axis=1)
        
    # Log transform data
    if log_transform:
        log2 = np.log2(data2)
        log2 = pd.concat([log2,data['Group']], axis=1)
        data = log2
        # print("NOTE we are Log2 ing this data, is that what you want?")
        # embed()
    else:
        data2['Group'] = data['Group']
        data = data2
    # embed()
    # split case and controls
    control = data[data['Group']==control][data.columns[:-1]]
    case = data[data['Group']==case][data.columns[:-1]]

    # drop 0 counts:
    control_counts = control.count()
    nas = control_counts[(control_counts!=control_counts.max())&(control_counts>0)]
    na_data = control[nas.index].dropna() 
    na_zscores = pd.DataFrame(stats.zscore(na_data, ddof=1), index=na_data.index, columns=na_data.columns)
    

    # calculate z-scores
    control_z = pd.DataFrame(stats.zscore(control, ddof=1), index=control.index, columns=control.columns)
    control_mean = control.mean()
    control_std = control.std()
    
    control_z[nas.index]=na_zscores[nas.index] 
    
    
    
    case_z = (case-control_mean)/control_std
    # embed()
    # only return cols with data:
    show=[]
    for c in data.columns[:-1]:
        if str(data[c].mean())!='nan':
            show.append(c)
    # embed()
            
    # control_z[show].to_csv('GWI_Control_zscore.csv')
    # case_z[show].to_csv('GWI_GWS_zscore.csv')

    
    return {'case':case_z[show], 'control':control_z[show]}
    
    
def make_VR_spreadsheet(infile, type, case, control, method, study, db_file='', table1_file='', FDR_threshold=0.05):
    data = pd.read_csv(infile,index_col=0)
    zscores = get_zscores(data, case, control)

    # data for CASES:
    zscores = zscores[type]
    # data for CONTROLS:
    
    make_corr_spreadsheet(zscores, type, method, infile, study, table1_file=table1_file, FDR_threshold=FDR_threshold)

def get_correlations(data, method='pearson'):
    ''' EXPECTS Z-scored data '''
    # embed()
    corr = data.corr(method=method)
    # embed()
    out = []
    idx=0
    anticorr=[]
    for i in corr.index:
        for c in corr.columns:
            name = []
            name.append(i)
            name.append(c)
            name.sort()
            name_out='|'.join(name)
            if method=='pearson':
                r, p = stats.pearsonr(data[c], data[i])
                
            elif method=='spearman':
                r, p = stats.spearmanr(data[c], data[i])
            
            # out.append({'from':i,'to':c, 'corr':data.loc[i,c]})
            if r==-1:
                print('perfect anticorrelation?!', c, i)
                anticorr.append((c,i))
            out.append({'name':name_out, 'r':r, 'p':p})
        idx+=1
    out = pd.DataFrame(out)
    if len(anticorr)>1:
        print('you have %i anticorrelated pairs'%len(anticorr))
        print(anticorr)
    # print('test')
    # embed()
    

    print('list len:', len(out))
    out = out.drop_duplicates()
    print('drop dupes:', len(out))
    
    # NOTE drop dupes doesnt quite work for spearman - numerical imprecision seems to cause some dupes to remain:
    tmp = out['name'].drop_duplicates()
    out = out.ix[tmp.index]
    print('drop dupes by name:', len(out))
    # embed()
    out['from']=out['name'].str.split('|',expand=True)[0]
    out['to']=out['name'].str.split('|',expand=True)[1]
    out = out[out['r']<1]
    out = out[out['from']!=out['to']]
    # embed()
    print('drop identity:', len(out))
    out['abs(r)'] = abs(out['r'])
    del out['name']

    # FDR:
    from statsmodels.stats.multitest import fdrcorrection_twostage
    res = fdrcorrection_twostage(out['p'], alpha=0.05, method='bh')
    rej = res[0]
    pval_corr = res[1]

    # Qvalue:
    import qvalue
    qv = qvalue.estimate(out['p']) # Storey's q value

    out['FDR'] = pval_corr
    out['Q value'] = qv
    
    # embed()

    return out
    
    
    
    
def map_to_mdtb(corr, db_file):
    
    alias = pd.read_excel(db_file,'Alias')
    mrm = pd.read_excel(db_file,'MRM')
    chemical = pd.read_excel(db_file,'Chemical', index_col=0)
    pathway = pd.read_excel(db_file,'Pathway')

    data = pd.merge(alias, mrm, left_on='MRM Index', right_on='MRM Index')
    # embed()
    data = pd.merge(data, chemical, left_on='Chemical Index', right_on='Chemical Index')
    data = pd.merge(data, pathway, left_on='Pathway Index', right_on='Pathway Index')
    
    # embed()
    
    # below added to remove dupe aliases 8/6/21
    data = data[['MRM Name_x','Chemical Index', 'Pathway Name', 'Exo/Endo']]
    data = data.groupby('MRM Name_x').first().reset_index()
    # end addition

    out = pd.merge(corr, data[['MRM Name_x','Chemical Index', 'Pathway Name','Exo/Endo']], left_on='from', right_on='MRM Name_x', how='left')
    out = pd.merge(out, data[['MRM Name_x','Chemical Index', 'Pathway Name','Exo/Endo']], left_on='to', right_on='MRM Name_x', how='left')
    
    # embed()
    
    in_pwy = out[out['Pathway Name_x']==out['Pathway Name_y']] 
    out_pwy = out[out['Pathway Name_x']!=out['Pathway Name_y']] 
    
    out.loc[in_pwy.index, 'Pathway']='IN'
    out.loc[out_pwy.index, 'Pathway']='OUT'
         
    del out['MRM Name_x_x']
    del out['MRM Name_x_y']
    
    
    out = out.rename(columns={'Chemical Index_x':'FROM CINDEX',
        'Chemical Index_y':'TO CINDEX',
        'Pathway Name_x':'Pathway1',
        'Pathway Name_y':'Pathway2',
        'Exo/Endo_x':'Endo/Exo1',
        'Exo/Endo_y':'Endo/Exo2',
    })
    
    return out
    
    
def map_to_cytoscape(data, cyto_xy, corr_type):
    # CID_MAP.to_excel('GWI_%s_%s_corr_list_ZSCORE_CHEM_INDEX.xlsx'%(type,corr_type))

    # embed()
    # make "VR" sheet
    cytoscape = pd.read_excel(cyto_xy)

    out2 = pd.merge(data, cytoscape, left_on='FROM CINDEX', right_on='Chemical Index', how='left')
    
    out2 = pd.merge(out2, cytoscape, left_on='TO CINDEX', right_on='Chemical Index', how='left', suffixes=('1','2'))
    
    no_map1=out2[out2['Chemical Index1'].isnull()][['from','FROM CINDEX']].drop_duplicates()
    no_map1 = no_map1.rename(columns={'from':'name','FROM CINDEX':'CINDEX'})
    no_map2=out2[out2['Chemical Index2'].isnull()][['to','TO CINDEX']].drop_duplicates()
    no_map2 = no_map2.rename(columns={'to':'name','TO CINDEX':'CINDEX'})
    no_map = no_map1.append(no_map2)
    no_map = no_map.drop_duplicates()
    # check how many pairwise connections there are:
    tmp2=out2.groupby(['from','to']).count().reset_index()[['from','to','x1']]
    tmp2 = tmp2.rename(columns={'x1':'Duplicate Pairs'})
    out2 = pd.merge(out2, tmp2, left_on=['from','to'], right_on=['from','to'])
    # embed()
    

    del out2['from']
    del out2['to']
    del out2['FROM CINDEX']
    del out2['TO CINDEX']

    out2 = out2.rename(columns={'p':'%s p'%corr_type, 'r':'%s r'%corr_type})
    # out2['Pair ID'] = range(1,len(out2)+1)

    out2 = out2.drop_duplicates() # NOTE several dupes emerge because the map has duplicate nodes
    
    out2 = out2.dropna() # NOTE this removes all correlations that ARENT on the maps
    
    return {'cyto':out2,'no_map':no_map}
    