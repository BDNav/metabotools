
import pandas as pd
import numpy as np
import pandas as _pandas
from collections import OrderedDict as _OrderedDict
import re as _re
import itertools as _itertools
import scipy.stats as _stats
import numpy as _numpy
import matplotlib.pyplot as _plot
import datetime as _datetime
import os as _os
import pandas as pd
from IPython import embed

from .metabotools import db_merge_full, db_merge_mrm, db_merge_name, db_merge_path, db_merge_chem, db_read_excel    


def main():
    # dir = 'test_data/Plasma_M/'
    # data_file = dir+'Plasma_M_suicide_Sho_rap_GDF+FGF_20190926_N=93.csv'
    # exp = 'Depress'
    # control = 'Control'
    # study_name = 'Plasma_M_Suicide_Study_102020'
    # cutoff=0
    # vip_file = dir+'plsda_vip.csv'
    # mda_file = dir+'randomforests_sigfeatures.csv'
    
    
    dir = 'test_data/Plasma_F/'
    data_file = dir+'Plasma_F_suicide_Sho_Rap_FGF+GDF_20190929_N=99.csv'
    exp = 'Depress'
    control = 'Control'
    study_name = 'Plasma_F_Suicide_Study_102020'
    cutoff=0
    vip_file = dir+'plsda_vip.csv'
    mda_file = dir+'randomforests_sigfeatures.csv'
    
    
    db_version = '6.2_protein'
    db_file = r'D:\Dropbox (Personal)\naviauxlab_informatics\mtdb\MasterDatabase_v%s.xlsx'%db_version


    md_num=1
    at='v'
    cutoff=0
    table1_out = dir+'%s_table_ONE_mtb%s_%s_%0.1f_MWQ_FILTERED.xlsx'%(study_name, db_version,at,cutoff)
    print('\t', table1_out)
    buildTableOne(data_file, exp, control, vip_file, md_num, 'v', cutoff, db_file, table1_out, mda_file=mda_file)
    
    embed()
    

class TableOne(object):
    
    def __init__(self, metabolite_database, data, metric_column, change_column, total_MRMs=1000.0):
        self.database = metabolite_database
        self.data = data
        self.metric_column = metric_column
        self.change_column = change_column
        self.dataset = db_merge_mrm(self.data, self.database, how='left')       
        self.total_MRMs = total_MRMs * 1.0
        self.table = self._construct_table()
            
    headers = _OrderedDict([('measure', u'Measured Metabolites in the Pathway (N)'),
               ('exp_prop', u'Expected Pathway Proportion (P = N/1006)'),
               ('exp_hits', u'Expected Hits in Sample of 61 (P * 61)'),
               ('obs', u'Observed Hits in the Top 61 Metabolites'),
               ('enrich1', u'Fold Enrichment (Obs/Exp)'),
               ('enrich2', u'Hypergeometric p Value'),
               # ('enrich3', u'Hypergeometric p Value EXACT'),
               ('impact', u'Impact (Sum X Score)'),
               ('frac_impact', u'Fraction of Impact (Z Score) Explained (% of 142.4701)'),
               ('inc', u'Increased'),
               ('dec', u'Decreased')])    
   
    def _construct_table(self):            
        grouped = self.dataset.groupby(['Pathway Index', 'Pathway Name'])
        
        tableOne = _pandas.DataFrame({self.headers['impact']: grouped.sum()[self.metric_column]})
        total_score = tableOne.sum()[0]            
        self.headers['frac_impact'] = u'Fraction of Impact ({}) Explained (% of {})'.format(self.metric_column, total_score)
        tableOne[self.headers['frac_impact']] = tableOne[self.headers['impact']].div(total_score)
       
        total_obs = grouped.size().sum()
        self.headers['obs'] = u'Observed Hits in the Top {} Metabolites'.format(total_obs)
        tableOne[self.headers['obs']] = grouped.size()
        tableOne[self.headers['dec']] = grouped.apply(lambda x: len(x[x[self.change_column] < 0]))
        tableOne[self.headers['inc']] = grouped.apply(lambda x: len(x[x[self.change_column] > 0]))

        c_groups = self.database.groupby(['Pathway Index', 'Pathway Name'])
        c_groups_size = c_groups.size()
        c_groups_size.name = self.headers['measure']
 
        tableOne = tableOne.join(c_groups_size, how='left')
        self.headers['exp_prop'] = u'Expected Pathway Proportion (P = N/{})'.format(int(self.total_MRMs))
        tableOne[self.headers['exp_prop']] = tableOne[self.headers['measure']].div(self.total_MRMs)
        self.headers['exp_hits'] = u'Expected Hits in Sample of {0} (P * {0})'.format(total_obs)
        tableOne[self.headers['exp_hits']] = tableOne[self.headers['exp_prop']].mul(total_obs)
        tableOne[self.headers['enrich1']] = tableOne[self.headers['obs']] / tableOne[self.headers['exp_hits']]
        
        # tableOne[self.headers['enrich2']] = tableOne.apply(lambda x: 1 - _stats.hypergeom.cdf(x[self.headers['obs']] - 1,
                                                                                             # int(self.total_MRMs),
                                                                                             # x[self.headers['measure']],
                                                                                             # total_obs), axis=1)
        # JON replaced above for CDF with sf, they are equivalent, sf is just 1-cdf
        tableOne[self.headers['enrich2']] = tableOne.apply(lambda x: _stats.hypergeom.sf(x[self.headers['obs']] - 1,
                                                                                             int(self.total_MRMs),
                                                                                             x[self.headers['measure']],
                                                                                             total_obs), axis=1)
        # JON ADDED THIS BELOW TO CALC EXACT P VALUE INSTEAD OF CDF:
        # see nice explanation here: https://blog.alexlenail.me/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
        # and a calculator that does both here: https://www.geneprof.org/GeneProf/tools/hypergeometric.jsp
        # change below to enrich3 and uncomment enrich3 line above to show both pvalues
        # tableOne[self.headers['enrich3']] = tableOne.apply(lambda x: _stats.hypergeom.pmf(x[self.headers['obs']],
                                                                                             # int(self.total_MRMs),
                                                                                             # x[self.headers['measure']],
                                                                                             # total_obs), axis=1)
        
        
                                                                                             
        return tableOne[list(self.headers.values())]
        '''
        res = tableOne[list(self.headers.values())]
        
        
        # scipy.stats.fisher_exact
        for i in res.index:
            x = res.loc[i,'Observed Hits in the Top {} Metabolites'.format(total_obs)]
            y = res.loc[i,'Measured Metabolites in the Pathway (N)']
            pro, p = _stats.fisher_exact([[x,y],[y,self.total_MRMs]])
            # print res.loc[i,'Pathway Name'],
            print p
            # if y==70:
                
            
            res.loc[i,'Fisher Exact p Value']=p
            
        
        return res
        '''
        
        
def buildTableOne(data_file, exp, control, vip_file, md_num, at, metric_cutoff, db_file, out_file, Q_filter_val=1, U_filter_val=1, mda_filter=None, mda_file='', replace_chars=True):
    '''
    replace_chars: replace (,:,) chars from metaboanalyst - these characters used to be removed by metaboanalyst, no longer the case for new versions
    '''
    analysis_type = {'z': 'Z Abs', 'v': 'VIP Score', 'm': 'Random Forest MDA (5000 trees)'}    
    db_dict = {'a': 'Alias',
           'c': 'Chemical',
           'm': 'MRM',
           'p': 'Pathway'}

    db_parts = db_read_excel(db_file, db_dict)
    db_parts['a'] = db_parts['a'].groupby('MRM Name').first() # Jon added 7/16/21 to remove problem of having duplicate aliases
    db = db_merge_full(db_parts['m'], db_parts['c'], db_parts['p'])
    # alias = pd.read_excel(db_file,'Alias')
    
    
    d = pd.read_csv(data_file)
    num_targeted = len(d.columns)-2
    d = d.dropna(axis=1, how='all').set_index(d.columns[0])
    
    
    g = pd.DataFrame(d[d.columns[0]])
    if g.columns[0]!='Group':
        print('Group must be second column')
    
    
    
    # check for any AUCS<0
    # embed()
    d2 = d[d.columns[md_num:]]<0
    if d2.fillna(0).sum().sum()<0:
        print('error, you have an AUC value < 0')
        d3=d2.fillna(0).sum()
        print('d3[d3<0]')
        embed()
    
    col = d.columns.tolist()
    m = d[col[0:md_num]]
    # data2 = d[col[md_num:]]
    # for c in data2.columns:
        # tmp = pd.to_numeric(data2[c])
        
    
    # embed()
    d = np.log2(d[col[md_num:]])
    col = pd.merge(pd.DataFrame(d.columns.tolist(), columns=['MRM Name']), db_parts['a'], on='MRM Name')
    # debug find missing cols:
    tmp = pd.merge(pd.DataFrame(d.columns.tolist(), columns=['MRM Name']), db_parts['a'], on='MRM Name', how='left')
    tmp = tmp.fillna('-')    
    missing = tmp[tmp['MRM Index']=='-']

    # dupes = pd.merge(tmp,db_parts['m'], left_on='MRM Index', right_on='MRM Index')
    # dupes = pd.merge(dupes, db_parts['c'], left_on='Chemical Index', right_on='Chemical Index')

    if len(missing)>0: # removed 12/3/21 caused errors on py3+
        print('ERROR: some measured mets dont have an alias')
        print(missing)
        # embed()
        # col = pd.merge(pd.DataFrame(d.columns.tolist(), columns=['MRM Name']), db_parts['a'], on='MRM Name')
        # col = col.fillna('-')
        # col = col[col['MRM Index']!='-']
        # d = d[col['MRM Name']]
    
    if len(col['MRM Index'])>len(d.columns):
        print('error MRM index merge ~= # columns: dupes in merge')
        tmp=col.groupby('MRM Name').count()
        print(tmp[tmp>1].dropna())
        col = col.sort_values(['MRM Name','MRM Index'])
        # START tmp addition
        tmp2 = col.groupby('MRM Name').last() # prioritize the UNX (U < T) markers
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.rename(columns={'index':'MRM Name'})
        col = tmp2
        embed()
        print('temp: addition for tox')
    
    col.index=col['MRM Name']
    
    # 1/22/2021 change to explicitly change header name not just rely on order
    d_rename = col['MRM Index']    
    d = d.rename(columns=d_rename.to_dict())
    
    # d.columns = col['MRM Index']
    # End 1/22/2021 change
    
    # g = pd.read_csv(group_file).set_index('Sample')
    # g = pd.read_csv(group_file).set_index('Sample')
    # embed()
    # db = db[db['MRM Index'].isin(d.columns.tolist())]
    db = db[db['MRM Index'].isin(d_rename.tolist())]
    # embed()
    
    z_score = pd.DataFrame(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()).reset_index().values.copy(), columns=['MRM Index', 'Z Score'])
    z_score['Z Abs'] = np.abs(z_score['Z Score'])
    

    geomean_ratio = pd.DataFrame(np.exp2(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean())).reset_index().values.copy(), columns=['MRM Index', 'Geo. Mean Ratio']) # old taylor code
    
    linear_ratio = pd.DataFrame(np.exp2(d.loc[g.index.tolist()][g.Group == exp]).mean().div(np.exp2(d.loc[g.index.tolist()][g.Group == control]).mean()).reset_index().values.copy(), columns=['MRM Index', 'Linear Ratio']) # jon update 9/15/2020
    
    # embed()
    
    
    samplesC = d.loc[g.index.tolist()][g.Group == exp]
    controlsC = d.loc[g.index.tolist()][g.Group == control]
    len_controls = len(controlsC)
    len_samples = len(samplesC)
    
    
    from scipy import stats 
    
    mann_whitney = {}
    mann_whitney2 = {}
    ttest = {}
    
    for c in controlsC.columns:
        controlC = controlsC[c].dropna()
        sampleC = samplesC[c].dropna()
        # if len(controlC)>len_controls*0.7 and len(sampleC)>len_samples*0.7:
        
        rej, pval = stats.ttest_ind(controlC, sampleC)
        ttest[c]=pval
        
        # Mann Whitney https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
        stat, p_mannw2 = stats.mannwhitneyu(controlC,sampleC, use_continuity=False, alternative='two-sided') # NOTE: these parameters will make the MWpval=KWpval (for two samples)
        stat, p_mannw = stats.mannwhitneyu(controlC,sampleC) # NOTE: these parameters will make the MWpval=KWpval (for two samples)
        mann_whitney[c]=p_mannw
        mann_whitney2[c]=p_mannw2
    mann_whitney = pd.Series(mann_whitney,name='Mann-Whitney p (1-sided)')
    mann_whitney2 = pd.DataFrame(pd.Series(mann_whitney2,name='Mann-Whitney U test p'))
           
    pvals = pd.DataFrame(pd.Series(ttest).dropna(),columns=['P Value'])
    
    ''' calculating storey's qvalue  '''
    # uses package: https://github.com/nfusi/qvalue'''
    from . import qvalue
    # embed()
    try:
        qv = qvalue.estimate(pvals.values) # Storey's q value
    except:
        print('Q value error - has occured when duplicate metabolites were present - at the unique ID stage')
        embed()
    qv = pd.DataFrame(qv, columns=['q-Value'], index=pvals.index)
    
    new_stats = pd.merge(mann_whitney2, qv, left_index=True, right_index=True)
    
    # embed()
    if vip_file:
        v = _pandas.read_csv(vip_file)
        
        v = v[[v.columns[0], v.columns[1]]]
        v.columns = ['MRM Name_x', 'VIP Score']
        # embed()


        if replace_chars==True:
            col['VIP2MRM']=col['MRM Name'].str.replace('(','').str.replace(')','').str.replace(':','').str.replace('[','').str.replace(']','')
        else:
            col['VIP2MRM']=col['MRM Name']
        
        v2 = pd.merge(v,col, left_on='MRM Name_x', right_on='VIP2MRM', how='right')
        
        v = pd.merge(v,col, left_on='MRM Name_x', right_on='VIP2MRM')
        v = v.sort_values('MRM Name')
        # v = pd.concat([v,col], axis=1)
        if len(v)!=len(col):
            print('merging error on vip scores')
            print("NOT IN MTDB:")
            print(v2[v2['MRM Index'].isna()])
            print("NOT IN VIP LIST:")
            print(v2[v2['VIP Score'].isna()])
            embed()
        
        data = pd.merge(pd.merge(v, z_score, on='MRM Index'), linear_ratio, on='MRM Index').infer_objects()[['MRM Index', 'VIP Score', 'Z Score', 'Z Abs', 'Linear Ratio','MRM Name']] # added MRM Name on 1/5/2021
        
        # print('dropping NA - this was added for Tox T1')
        
        data = data.dropna()
    else:
        data = pd.merge(z_score, linear_ratio, on='MRM Index').infer_objects()[['MRM Index', 'Z Score', 'Z Abs', 'Linear Ratio']]
    
    
    if mda_file!='':
        mda = _pandas.read_csv(mda_file)
        mda = mda[[mda.columns[0], mda.columns[1]]]
        mda.columns = ['MRM Name_x', 'MDA']
        
        if replace_chars==True:
            col['VIP2MRM']=col['MRM Name'].str.replace('(','').str.replace(')','').str.replace(':','').str.replace('[','').str.replace(']','')
        else:
            col['VIP2MRM']=col['MRM Name']
        
        
        mda2 = pd.merge(mda,col, left_on='MRM Name_x', right_on='VIP2MRM', how='right')
        
        mda = pd.merge(mda,col, left_on='MRM Name_x', right_on='VIP2MRM')
        mda = mda.sort_values('MRM Name')
        # v = pd.concat([v,col], axis=1)
        if len(mda)!=len(col):
            print('merging error on vip scores')
            print("NOT IN MTDB:")
            print(mda2[mda2['MRM Index'].isna()])
            print("NOT IN VIP LIST:")
            # print mda2[mda2['VIP Score'].isna()]
            embed()
            
        data = pd.merge(data, mda[['MRM Index','MDA']], left_on='MRM Index', right_on='MRM Index')
        data = data.rename(columns={'MDA':'Random Forest MDA (5000 trees)'})
    
    data = pd.merge(data, geomean_ratio, left_on='MRM Index', right_on='MRM Index')
    data = pd.merge(data, new_stats, left_on='MRM Index', right_index=True)
    
    
    
    '''
            
        data = pd.merge(pd.merge(mda, z_score, on='MRM Index'), linear_ratio, on='MRM Index').infer_objects()[['MRM Index', 'MDA', 'Z Score', 'Z Abs', 'Linear Ratio']]
        
        
        
    else:
        
            data = pd.merge(data, geomean_ratio, left_on='MRM Index', right_on='MRM Index')
            
    '''
    # embed()
    if at == 'z':
        trunc_data = data[data['Z Abs'] > metric_cutoff].sort_values(by='Z Abs', ascending=False)
    elif at == 'v':
        trunc_data = data[data['VIP Score'] > metric_cutoff].sort_values(by='VIP Score', ascending=False)
    elif at=='m':
        trunc_data = data[data['Random Forest MDA (5000 trees)'] > 0].sort_values(by='Random Forest MDA (5000 trees)', ascending=False)
        
    trunc_data2 = trunc_data
    
    trunc_data2 = trunc_data2[trunc_data2['q-Value'] <=Q_filter_val]
    q_filter = len(trunc_data2)
    trunc_data2 = trunc_data2[trunc_data2['Mann-Whitney U test p'] <=U_filter_val]
    u_filter = len(trunc_data2)
    
    if mda_filter!=None:
        trunc_data2 = trunc_data2[trunc_data2['Random Forest MDA (5000 trees)'] >mda_filter]
    
    
    
    
    mwq_filter=trunc_data2 # to be saved as a new tab
    
    if len(trunc_data2)==0:
        print('your filters are too strict, there is no data left to analyze')
        embed()
    
    t1 = TableOne(db, trunc_data2, metric_column=analysis_type[at], change_column='Z Score', total_MRMs=len(db))

    t1 = t1.table
    t1 = t1.sort_values('Impact (Sum X Score)', ascending=False)
    t1 = t1.rename(columns={'Impact (Sum X Score)':'Impact (Sum %s)'%analysis_type[at]})
    t1 = t1.reset_index()
    
    cols = t1.columns
    t1 = t1.append(pd.DataFrame([{'Increased':t1.Increased.sum(), 'Decreased':t1.Decreased.sum()}]))
    t1 = t1[cols]
    
    
    if mda_file!='':
        trunc_data2 = trunc_data2[trunc_data2['Random Forest MDA (5000 trees)'] > 0].sort_values(by='Random Forest MDA (5000 trees)', ascending=False)
       
    
        t1_mda = TableOne(db, trunc_data2, metric_column=analysis_type['m'], change_column='Z Score', total_MRMs=len(db))
        t1_mda = t1_mda.table
        t1_mda = t1_mda.sort_values('Impact (Sum X Score)', ascending=False)
        t1_mda = t1_mda.rename(columns={'Impact (Sum X Score)':'Impact (Sum %s)'%'MDA'})
        t1_mda = t1_mda.reset_index()
        cols = t1_mda.columns
        t1_mda = t1_mda.append(pd.DataFrame([{'Increased':t1_mda.Increased.sum(), 'Decreased':t1_mda.Decreased.sum()}]))
        t1_mda = t1_mda[cols]
    
    
    
    metadata = [
        _datetime.datetime.now(),
        _os.path.abspath(db_file), 
        _os.path.abspath(data_file),
        metric_cutoff,
        len(g[g==control].dropna()),
        len(g[g==exp].dropna()),
        num_targeted,
        len(d.columns),
        exp,
        out_file.split('/')[-1].split('_table_ONE')[0],
        Q_filter_val,
        U_filter_val,
        str(mda_filter),
        ]
                
    meta_index = [
        'Date',
        'Database File',
        'Metabolite File',
        '%s Cutoff'%analysis_type[at],
        'N (%s)'%control,
        'N (%s)'%exp,
        'Metabolites Targeted',
        'Metabolites Measured',
        'Group',
        'Study',
        'Q-value <=',
        'Mann-Whitney U value <=',
        'MDA filter >='
        ]
    
            
    metadata = _pandas.DataFrame({'':meta_index, 'TableOne Analysis': metadata})
    if vip_file:
        headers = ['MRM Name', 'Pathway Name', 'VIP Score', 'Z Score', 'Z Abs', 'Linear Ratio', 'Geo. Mean Ratio','Mann-Whitney U test p','q-Value']
    else:
        headers = ['MRM Name', 'Pathway Name', 'Z Score', 'Z Abs', 'Linear Ratio', 'Geo. Mean Ratio','Mann-Whitney U test p','q-Value']
        
    if mda_file!='':
        headers +=['Random Forest MDA (5000 trees)']
        
    tmp=data.groupby('MRM Index').count()
    
    tmp.sort_values('VIP Score')
    dupes = tmp[tmp['VIP Score']>1]
    if len(dupes)>0:
        print('error duplicates',dupes)
        # one option is to do this: data.groupby('MRM Index').first().reset_index()
        # but it seems to lose some data
        data = data.groupby('MRM Index').first().reset_index()
        
    
    data = db_merge_mrm(data, db, how='left')
    data = data.rename(columns={'MRM Name_x':'MRM Name'})
    data = data[headers].sort_values('VIP Score', ascending=False)
    
    mwq_filter = db_merge_mrm(mwq_filter, db, how='left')
    mwq_filter = mwq_filter.rename(columns={'MRM Name_x':'MRM Name'})
    mwq_filter = mwq_filter[headers].sort_values('VIP Score', ascending=False)
        
    writer = _pandas.ExcelWriter(out_file)
    metadata.to_excel(writer, 'Header', index=False)
    data.to_excel(writer, 'Metabolites', index=False)
    mwq_filter.to_excel(writer, 'Metabolites filtered', index=False)
    t1.to_excel(writer, 'Table 1', index=False)
    if mda_file!='':
        t1_mda.to_excel(writer, 'Table 1 MDA', index=False)
    writer.save()       
    
    
if __name__=='__main__':
    main()








