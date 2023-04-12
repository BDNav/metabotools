import pandas as pd
from IPython import embed
import numpy as np

#script adapted from MRM\Bob 12-7-16\analyze_data.py

def main():
    
    outfile='Plasma_M_TableS1.xlsx'
    study_name = 'test'
    data_dir = r'D:\Dropbox (Personal)\naviauxlab_informatics\Suicide_Study\2_26_2021\Plasma_M'
    data_file = data_dir + r'\Plasma_M_suicide_Sho_rap_GDF+FGF_20190926_N=93.csv'
    vip_file = data_dir + r'\plsda_vip.csv'
    mda_file = data_dir + r'\randomforests_sigfeatures.csv'
    
    control = 'Control'
    exp = 'Depress'
    table1_mets = data_dir + r'\Plasma_M_Suicide_Study_112020_table_ONE_mtb6.2_protein_MWQ_FILTERED.xlsx'
    table1_mets = pd.read_excel(table1_mets,'Metabolites')
    
    db_version = '6.2_protein'
    db_file = r'D:\Dropbox (Personal)\naviauxlab_informatics\mtdb\MasterDatabase_v%s.xlsx'%db_version

    
    new_build_table_s1(outfile, study_name, data_file, control, exp, table1_mets, db_file, VIP_cutoff=0.9, mda_file=mda_file, run_het=False,       knn_vip_threshold=1, add_RF_AUCs=True, add_logistic_reg_AUCs=True)
    embed()
    
def new_build_table_s1(outfile, study_name, data_file, control, exp, table1_mets, db_file, VIP_cutoff='', mda_file='', run_het=False,       knn_vip_threshold=1, add_RF_AUCs=False, add_logistic_reg_AUCs=False, overwrite=1, replace_chars=True):
    '''
    overwrite = 1 will re-run RF and logred (which take a long time) if they havent already been generated
    
    '''
    import os
    
    out = {}
    data = pd.read_csv(data_file)

    data = data.transpose().dropna().transpose()
    
    
    num_mets_targeted = len(data.columns)-2 # assuming first two are "Sample" and "Group"
    
    db = merge_db(db_file)
    
    ''' ORGANIZE DATA '''
    
    controls=data[data['Group']==control]
    samples=data[data['Group']==exp]

    
    controls_linear = controls[controls.columns[2:]]
    samples_linear = samples[samples.columns[2:]]
    for c in controls_linear.columns:
        controls_linear.loc[:,c] = pd.to_numeric(controls_linear.loc[:,c])
    for c in samples_linear.columns:
        samples_linear.loc[:,c] = pd.to_numeric(samples_linear.loc[:,c])
    controlsC = controls_linear.apply(np.log2).replace([np.inf,-np.inf],np.nan)
    samplesC = samples_linear.apply(np.log2).replace([np.inf,-np.inf],np.nan)
    
    met_counts = data.count()
    to_keep = met_counts[met_counts>0]
    
    stats = calc_stats(samples, samples_linear, samplesC, controls, controls_linear, controlsC)
    
    if run_het==True:
        het = calc_heteroscasdicity()
        
    sem = calc_sem_etc(controls_linear, controlsC, samples_linear, samplesC)
    
    other_stats = calc_other_stats(controls, controls_linear, controlsC, samples, samples_linear, samplesC)
    
    outliers = get_outliers(controlsC, samplesC)
    
    # get zscores:
    data = pd.read_csv(data_file, index_col=0)
    zscores = get_zscores(data, exp, control)
    zscores = zscores['case'].append(zscores['control'])
    zscores['Group'] = data['Group']
    zscores = zscores[['Group'] + zscores.columns[:-1].tolist()]
    
    
    corrs = calc_corr(data, control, exp)
    
    
    #MERGE ALL new stats:
    out = pd.merge(stats, sem, left_index=True, right_index=True)
    num_mets_measured = len(out)
    # embed()
    
    out = pd.merge(out, other_stats, left_index=True, right_index=True)
    out = pd.merge(out, outliers, left_index=True, right_index=True)
    
    out = pd.merge(out, corrs, left_index=True, right_index=True)
    if len(out)!=num_mets_measured:
        print('error merging - we lost a metabolite 1')
        embed()
        
    if mda_file!='':
        kmeans_clusters = [5, 10]
        mda_res = calc_knn_kmeans(zscores, mda_file, control, exp, kmeans_clusters, replace_chars) # TODO
        out = pd.merge(out, mda_res, left_index=True, right_index=True)
        # embed()
    
    if add_RF_AUCs==True:
        filename = 'tmp_RF_AUC_'+study_name+'.csv'
        if os.path.isfile(filename) and overwrite==0:
            indiv_RF_res = pd.read_csv(filename, index_col=0)
        else:
            indiv_RF_res = calc_indiv_RF_res(controls_linear, samples_linear)
            indiv_RF_res.to_csv(filename)
        
        out = pd.merge(out, indiv_RF_res, left_index=True, right_index=True)
        if len(out)!=num_mets_measured:
            print('error merging - we lost a metabolite 2')
            embed()
            
    if add_logistic_reg_AUCs==True:
        filename = 'tmp_logreg_AUC_'+study_name+'.csv'
        if os.path.isfile(filename) and overwrite==0:
            indiv_logreg = pd.read_csv(filename, index_col=0)
        else:
            #indiv_logreg = calc_logreg(controlsC, samplesC)
            indiv_logreg = calc_logreg(controls_linear, samples_linear)
            indiv_logreg.to_csv(filename)
        # indiv_logreg = pd.read_csv('logreg.csv', index_col=0)
        out = pd.merge(out, indiv_logreg, left_index=True, right_index=True)
        if len(out)!=num_mets_measured:
            print('error merging - we lost a metabolite 3')
            embed()
        
    
    # RENAMES:
    out['Percent Difference (Case/Control -1)'] = out['Fold Change']-1
    out['Standard Error of the Difference (%)'] = out['Linear SEM']-1
    
    if 'Mann-Whitney U test p' in table1_mets.columns:
        del table1_mets['Mann-Whitney U test p']
        del table1_mets['q-Value']
    
    
    # MERGE WITH table1_mets
    out = pd.merge(table1_mets, out, left_on='MRM Name', right_index=True)
    # MERGE WITH MTDB
    del db['MRM Name']
    
    db = db[['MRM Alias','MRM Index','CAS Number', 'HMDB number']].drop_duplicates()
    
    out2 = pd.merge(out, db, left_on='MRM Name', right_on='MRM Alias') 
    out2 = out2.groupby('MRM Name').first()
    out2 = out2.reset_index()
    
    
    if len(out2)!=num_mets_measured:
        print('error merging - we lost a metabolite 4')
        embed()
    out = out2
    del out['MRM Alias']
    
    
    num_mets_measured = len(out)
    
    ''' BELOW IS FROM old TS1 - time to update '''
    '''
    1 merge with db
    2. merge all dataset and rename columns
    3. merge with table1 mets
    
    '''
    
    ''' merge it all together '''
    cols_to_include = [
        'MRM Name',
        'HMDB number',
        'CAS Number',
        'Pathway Name',
        'VIP Score',
        'Z Score',
        'Linear Z Score',
        'P Value',
        'Kruskal-Wallis p',
        # 'Mann-Whitney p (1-sided)', # don't included 1-sided per Bob Email 10/17/2018
        'Mann-Whitney U test p',
        'Welch\'s ttest p',
        'D\'Agostino & Pearson omnibus p',
        'D\'Agostino & Pearson omnibus p (Linear)',
        'Shapiro-Wilk normality p',
        'Shapiro-Wilk normality p (Linear)',
        'Control Outliers',
        'Control Low Outlier Threshold',
        'Control High Outlier Threshold',
        'Case Outliers',
        'Case Low Outlier Threshold',
        'Case High Outlier Threshold',
        'FDR',
        'q-Value',
        'Control Mean',
        'Control Median',
        'Control SD*',
        'Control SEM*',
        'Case Mean',
        'Case Median',
        'Case SD*',
        'Case SEM*',
        'Linear Control Mean (Geomean)',
        '2.50%','97.50%','Hi/Lo Ratio',
        'Linear Case Mean (Geomean)',
        'Fold Change',
        'Log2 SEM',
        'Linear SEM',
        'Percent Difference (Case/Control -1)',
        'Standard Error of the Difference (%)',
        'Spearman_r',
        'Spearman_p',
        'Pearson_r', 
        'Pearson_p',
        'Kendall_tau',
        'Kendall_p',
        ]
        
    if add_RF_AUCs==True:
        cols_to_include+=indiv_RF_res.columns.tolist()
        
    if add_logistic_reg_AUCs==True:
        cols_to_include+=indiv_logreg.columns.tolist()
        
    if run_het==True:
        cols_to_include += [
            'Control Heteroscedasticity (linear space Breusch-Pagan p)',
            'Control Heteroscedasticity (log space Breusch-Pagan p)',
            'Case Heteroscedasticity (linear space Breusch-Pagan p)',
            'Case Heteroscedasticity (log space Breusch-Pagan p)'
        ]   
    
    
    
    if mda_file!='':
        cols_to_include.append('RFA MDA (5000 trees)')
        for n_clusters in kmeans_clusters:
            cols_to_include.append('kNN Cluster (of %i)'%n_clusters)
            cols_to_include.append('MDA Rank in kNN Cluster (of %i)'%n_clusters)
            cols_to_include.append('Mean MDA/kNN Cluster (of %i)'%n_clusters)
            cols_to_include.append('K-Means Cluster (of %i)'%n_clusters)
            
    
    out = out[cols_to_include] 
    
    out.index = range(1,len(out)+1)

    import datetime as dt
    if VIP_cutoff=='':
        print('new code calls for listing VIP cutoff - see CFS2 run all')
        embed()
        
    # RENAME COLUNS PER BOB's request 9/24/2018    
    out = out.rename(columns={
        'VIP Score':'PLSDA VIP Score',
        'P Value':'ttest p',
        'Control Mean':'Log2 Control Mean',
        'Control Median':'Log2 Control Median',
        'Case Mean':'Log2 Case Mean',
        'Case Median':'Log2 Case Median',
        '2.50%':'2.50% Threshold#',
        '97.50%':'97.50% Threshold+'
    })
    # embed()
    if mda_file!='':
        out = out.rename(columns={
            'RFA MDA (5000 trees)':'Random Forest MDA (5000 trees)'
        })
        tmp = out['Random Forest MDA (5000 trees)'].sort_values(ascending=False)
        tmp = pd.DataFrame(tmp)
        tmp['Random Forest Rank'] = range(1,len(tmp)+1)
        del out['Random Forest MDA (5000 trees)']
        out = pd.merge(out, tmp, left_index=True, right_index=True)
        
    
    # renaming cols per Bob's request - 2-11-2019
    out = out.rename(columns={
        'Z Score':'Log Space Z-Score',
        'ttest p':'Student\'s ttest p',
        # 'kNN Cluster':'kNN Cluster (of 10)'
    })
    
    writer = pd.ExcelWriter(outfile)
    readme = {'Study':study_name,
    'Group':exp,
    'Database File':db_file,
    'Analysis date':dt.datetime.today().strftime("%m/%d/%Y"),
    'Input Filename':data_file,
    'VIP Score Cutoff':VIP_cutoff,
    'N (control=%s)'%control:len(controls),
    'N (case=%s)'%exp:len(samples),
    'Metabolites Measured':num_mets_measured,
    'Metabolites Targeted':num_mets_targeted,
    'Methods':'Outlier calculation was by the method of Pierce implemented in Python',
    'Methods ':'Log2 SEM calculation was by dividing the standard deviation of the original distribution by the square root of the sample size for Log2 transformed data',
    'Methods  ':'Linear SEM calculation was by dividing the standard deviation of the original distribution by the square root of the sample size for Linear data',
    'Notes':'All p values are for 2-sided tests unless otherwise specified',
    'knn VIP cutoff':knn_vip_threshold
    # 'ML iterations':ML_iter,
    }
    readme = pd.DataFrame(pd.Series(readme),columns=['Table S1'])
    readme.to_excel(writer,'README')
    
    
    
    zscores.to_excel(writer,'Zscores')
    out.to_excel(writer,'TableS1')
    
    if len(out)!=len(zscores.columns)-1:
        print('tables1 lenght != zscores - something got dropped...')
        embed()
        
    
    if mda_file!='':
        knn_mda_files, knn_mda_vip_files = code_to_package_knn_mda(mda_res, out, kmeans_clusters, outfile, knn_vip_threshold) # TODO don't forget to package this as well
        # embed()
        for n_clusters in kmeans_clusters:
            knn_mda_files[str(n_clusters)].to_excel(writer,'kNN (%s) Ranking by MDA'%str(n_clusters))
            knn_mda_vip_files[str(n_clusters)].to_excel(writer,'kNN (%s) Ranking by VIP>%0.1f'%(n_clusters,knn_vip_threshold))
    writer.save()
    
    
    
def merge_db(db_file):
    mrm = pd.read_excel(db_file,'MRM')
    alias = pd.read_excel(db_file,'Alias')
    alias = alias.rename(columns={'MRM Name':'MRM Alias'})
    chem = pd.read_excel(db_file,'Chemical')
    db = pd.merge(mrm,alias, left_on='MRM Index', right_on='MRM Index')[['MRM Name','MRM Index','MRM Alias','Chemical Index']]
    db = pd.merge(db,chem[['Chemical Index','CAS Number','HMDB number']], left_on='Chemical Index', right_on='Chemical Index', how='left')
    
    return db
    
    

def calc_corr(data, control, exp):
    ''' ADD CORRELATIONS '''
    from scipy import stats
    data_col=1
    
    data2 = {}
    for c in data.columns[data_col:]:
        try:
            data2[c] = pd.to_numeric(data[c])
        except:
            print('error', c)
            embed()
    data2 = pd.DataFrame(data2, index=data.index)
    data2 = np.log(data2)
    # embed()
    groups = data['Group']
    groups[groups==exp]=1
    groups[groups==control]=0
    from scipy.stats import spearmanr, pearsonr
    corr_out = {}
    
    for c in data2.columns:
      
        if str(data2[c].mean())!='nan':
            vals = data2[c].dropna()
            
            group = groups.loc[vals.index]
            
            # adding kendall tau correlation (per Bob request 11/5/2020
            ktau, kpval = stats.kendalltau(group.values,vals.values) # note I am arbitrarily cutting off the last item... may not be fair
            
            corr_out[c] = {
                'Spearman_r': spearmanr(group.values, vals.values)[0], 
                'Spearman_p': spearmanr(group.values, vals.values)[1],
                'Pearson_r': pearsonr(group.values, vals.values)[0], 
                'Pearson_p': pearsonr(group.values, vals.values)[1],
                'Kendall_tau': ktau, 
                'Kendall_p': kpval
            }
            
    corr_out = pd.DataFrame(corr_out).transpose()
    return corr_out
    

def code_to_package_knn_mda(knn_mda, out, kmeans_clusters, outfile, knn_vip_threshold):
    def package_mdas(knn_mda_in, n_clusters):
        # embed()
        knn_mda2 = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).sum()
        knn_mda2['Count per cluster'] = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']
        positive_mdas = knn_mda_in[knn_mda_in['Random Forest MDA (5000 trees)']>0]
        # embed()
        knn_mda2['Count of positive MDAs per cluster'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']
        knn_mda2['Mean Z-score'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).mean()['Log Space Z-Score']
        knn_mda2['SD'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).std()['Log Space Z-Score']
        # embed()
        knn_mda2['Sum of positive MDAs per cluster'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).sum()['Random Forest MDA (5000 trees)']
        knn_mda2['Fraction of positive']=knn_mda2['Sum of positive MDAs per cluster']/knn_mda2['Sum of positive MDAs per cluster'].sum()
        
        knn_mda2 = knn_mda2.sort_values('Fraction of positive', ascending=False)
        knn_mda2['Rank']=range(1,len(knn_mda2)+1)
        
        for i in range(1,11):
            cluster_mdas = positive_mdas[positive_mdas['kNN Cluster (of %i)'%n_clusters]==i]
            cluster_mdas = cluster_mdas.sort_values('Random Forest MDA (5000 trees)', ascending=False)[:20]
            # embed()
            count=1
            
            for ii in cluster_mdas.index:
                name = cluster_mdas.loc[ii,'MRM Name']
                knn_mda2.loc[i,count] = name
                count+=1
        knn_mda2 = knn_mda2[[u'Random Forest MDA (5000 trees)', u'Log Space Z-Score',
           u'Count per cluster', u'Count of positive MDAs per cluster',
           u'Sum of positive MDAs per cluster',
           u'Fraction of positive', u'Mean Z-score', u'SD', u'Rank']+list(range(1,21))]
        return knn_mda2
        
    def package_knn_vip(knn_mda_in, n_clusters):
        # embed()
        knn_mda2 = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).sum()
        knn_mda2['Count per cluster'] = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']

        positive_mdas = knn_mda_in[knn_mda_in['PLSDA VIP Score']>knn_vip_threshold]
        # embed()
        knn_mda2['Mean Z-score'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).mean()['Log Space Z-Score']
        knn_mda2['SD'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).std()['Log Space Z-Score']
        # positive_mdas = knn_mda_in
        knn_mda2['Count of VIP>%0.1f per cluster'%(knn_vip_threshold)] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']
        # embed()
        knn_mda2['Sum of VIP>%0.1f per cluster'%(knn_vip_threshold)] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).sum()['PLSDA VIP Score']

        knn_mda2['Fraction of positive']=knn_mda2['Sum of VIP>%0.1f per cluster'%(knn_vip_threshold)]/knn_mda2['Sum of VIP>%0.1f per cluster'%(knn_vip_threshold)].sum()

        knn_mda2 = knn_mda2.sort_values('Fraction of positive', ascending=False)
        knn_mda2['Rank']=range(1,len(knn_mda2)+1)
                
        for i in range(1,11):
            cluster_mdas = positive_mdas[positive_mdas['kNN Cluster (of %i)'%n_clusters]==i]
            cluster_mdas = cluster_mdas.sort_values('PLSDA VIP Score', ascending=False)[:20]
            # embed()
            count=1
            for ii in cluster_mdas.index:
                name = cluster_mdas.loc[ii,'MRM Name']
                knn_mda2.loc[i,count] = name
                count+=1
        del knn_mda2['Log Space Z-Score']
        
        num_clusters = len(knn_mda2.columns)-8
        
        knn_mda2 = knn_mda2[[u'PLSDA VIP Score', u'Count per cluster',
           u'Count of VIP>%0.1f per cluster'%(knn_vip_threshold), u'Sum of VIP>%0.1f per cluster'%(knn_vip_threshold),
           u'Fraction of positive', 'Mean Z-score','SD',u'Rank']+list(range(1,num_clusters+1))]
     
        return knn_mda2
        
    
    knn_mda_files = {}
    knn_mda_vip_files = {}
    for n_clusters in kmeans_clusters:
        knn_outfile = outfile.replace('TableS1','knnClusters_%i'%n_clusters)
        writer2 = pd.ExcelWriter(knn_outfile)
        
        knn_mda = out[['MRM Name','Random Forest MDA (5000 trees)','kNN Cluster (of %i)'%n_clusters,'Log Space Z-Score']]
        
        # addition 7/23/2020 to made a new table for same analysis but with vip>1
        knn_mda_vip = out[['MRM Name','kNN Cluster (of %i)'%n_clusters,'PLSDA VIP Score','Log Space Z-Score']]# filter for those with vip>1
        # knn_mda_vip = knn_mda_vip[knn_mda_vip['PLSDA VIP Score']>1]
        # embed()
        knn_mda_vip_packaged = package_knn_vip(knn_mda_vip, n_clusters)
        # end new analysis 7/23/2020
            
        knn_mda_packaged = package_mdas(knn_mda, n_clusters)
        knn_mda_files[str(n_clusters)] = knn_mda_packaged
        knn_mda_vip_files[str(n_clusters)] = knn_mda_vip_packaged
        
        out['kNN Cluster (of %i)'%n_clusters] = out['kNN Cluster (of %i)'%n_clusters].fillna(99)
        # embed()
            
        knn_groups = out['kNN Cluster (of %i)'%n_clusters].unique().tolist()
        knn_groups.sort()
        # embed()
        for knn_group in knn_groups:
            groupmets = out[out['kNN Cluster (of %i)'%n_clusters]==knn_group]
            groupmets = groupmets[['MRM Name','PLSDA VIP Score','Log Space Z-Score','Mann-Whitney U test p','FDR','Random Forest MDA (5000 trees)','Pathway Name','kNN Cluster (of %i)'%n_clusters]]
            groupmets.to_excel(writer2, str(int(knn_group)))
       
        writer2.save()
    return knn_mda_files, knn_mda_vip_files
    
def calc_knn_kmeans(zscores, mda_file, control, exp, kmeans_clusters, replace_chars=True):


    new_analysis = pd.DataFrame([], index=zscores.columns[1:])  
    new_analysis['metaboanalyst_name']=new_analysis.index
    
    if replace_chars==True:
        new_analysis['metaboanalyst_name'] = new_analysis['metaboanalyst_name'].str.replace(':','')
        new_analysis['metaboanalyst_name'] = new_analysis['metaboanalyst_name'].str.replace('(','')
        new_analysis['metaboanalyst_name'] = new_analysis['metaboanalyst_name'].str.replace(')','')
    
    
    
    ''' MEAN DECREASE IN ACCURACY '''
    # TOO SLOW NOW in PYTHON
    # USE results from metaboanalyst for now
    
    mda = pd.read_csv(mda_file, index_col=0)
    mda = mda.rename(columns={'MeanDecreaseAccuracy':'RFA MDA (5000 trees)'})
        
    new_analysis = pd.merge(new_analysis, mda, left_on='metaboanalyst_name', right_index=True)
    
    if len(new_analysis)!=len(mda):
        print('error merging metaboanalysis MDA file with our metabolites')
        # check new_analysis metaboanalyst name above
        embed()
    
    ''' KMEANS ADDITIONS '''
        
    from sklearn.cluster import KMeans
    
    X=zscores[zscores.columns[1:]].dropna().transpose() # drops any samples with NA values
    
    
    all_kmeans=[]
    for n_clusters in kmeans_clusters:
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
        y_kmeans = kmeans.predict(X)
        km_out = pd.Series(y_kmeans,index=X.index)
        km_out+=1
        km_out.name = 'K-Means Cluster (of %i)'%n_clusters
        
        all_kmeans.append(km_out)
    all_kmeans=pd.DataFrame(all_kmeans).transpose()
    
    
    
    ''' KNN Additions '''
    if mda_file!='':
        X = zscores[zscores.columns[1:]].dropna().transpose()
        
        for n_clusters in kmeans_clusters:
            from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
            from sklearn.cluster import SpectralClustering
            y = all_kmeans['K-Means Cluster (of %i)'%n_clusters]
            knn = KNeighborsClassifier(n_neighbors=n_clusters).fit(X,y)
            knn = pd.DataFrame(knn.predict(X), index=X.index)
            
            knn=knn.rename(columns={0:'kNN Cluster (of %i)'%n_clusters})
            
            new_analysis_tmp = pd.merge(new_analysis, knn, left_index=True, right_index=True, how='left')
            all_ranks = pd.DataFrame()
            tmp = new_analysis_tmp[['RFA MDA (5000 trees)','kNN Cluster (of %i)'%n_clusters]]
            
            for cluster in tmp['kNN Cluster (of %i)'%n_clusters].unique():
                tmp2 = tmp[tmp['kNN Cluster (of %i)'%n_clusters]==cluster]
                tmp2=tmp2.sort_values('RFA MDA (5000 trees)', ascending=False)
                tmp2['Rank']=range(1,len(tmp2)+1)
                all_ranks = pd.concat([all_ranks,tmp2])
            tmp2 = tmp.groupby('kNN Cluster (of %i)'%n_clusters).mean()
            all_ranks = all_ranks.rename(columns={'Rank':'MDA Rank in kNN Cluster (of %i)'%n_clusters})
            tmp2 = tmp2.rename(columns={'RFA MDA (5000 trees)':'Mean MDA/kNN Cluster (of %i)'%n_clusters})
            all_ranks = pd.merge(all_ranks,tmp2, left_on='kNN Cluster (of %i)'%n_clusters, right_index=True)
            
            del all_ranks['RFA MDA (5000 trees)']
            new_analysis = pd.merge(new_analysis, all_ranks, left_index=True, right_index=True, how='left')
    new_analysis = pd.merge(new_analysis, all_kmeans, left_index=True, right_index=True)
            
    return new_analysis
    
def calc_indiv_RF_res(controlsC, samplesC, num_trees=100):
    ''' ADDING RF 3/4/21 '''
    
    tmp_controlsC = controlsC.copy()
    tmp_controlsC['pheno'] = 0
    tmp_samplesC = samplesC.copy()
    tmp_samplesC['pheno'] = 1

    data = tmp_controlsC.append(tmp_samplesC)
    cnts = data.count()
    to_keep = cnts[cnts>0]
    data = data[to_keep.index]
    
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score
    from sklearn import metrics
    from sklearn.model_selection import cross_val_score
    from scipy import stats
    
    
    out = {}
    
    for c in data.columns[:-1]:
        combo = data[[c,'pheno']].dropna()
        case = combo[combo['pheno']==1][c]
        control = combo[combo['pheno']==0][c]
        
        res, pval = stats.ttest_ind(case, control)
        
        
        X = combo[[c]]
        y = combo[['pheno']]
        
        
        # X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.2, random_state=42)
        
        clf = RandomForestClassifier(max_depth=2, random_state=0, n_estimators=num_trees)
        result = clf.fit(X, y)
        # scores = cross_val_score(clf, X, y, cv=3, scoring='roc_auc', n_jobs=6)
        scores = cross_val_score(clf, X, y, cv=3, scoring='roc_auc')
        
        out[c]={
            'AUROC (Random Forest)':scores.mean(),
            'RF AUROC Lower Limit (95%CI)':scores.min(),
            'RF AUROC Upper Limit (95%CI)':scores.max(),
            # 'ttest-pval':pval
        }
    out = pd.DataFrame(out).transpose()
    
    return out
    
def calc_logreg(controlsC, samplesC):
    ''' ADDING log reg 2/27/21 '''
    
    tmp_controlsC = controlsC.copy()
    tmp_controlsC['pheno'] = 0
    tmp_samplesC = samplesC.copy()
    tmp_samplesC['pheno'] = 1

    logreg_data = tmp_controlsC.append(tmp_samplesC)
    logreg_cnts = logreg_data.count()
    to_keep = logreg_cnts[logreg_cnts>0]
    logreg_data = logreg_data[to_keep.index]
    
    import statsmodels.api as sm
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import cross_val_score
    logreg_out = {}
    random_states = [10, 42, 100, 33, 12]
    
    for c in logreg_data.columns[:-1]:
        aurocs = []
        combo = logreg_data[[c,'pheno']].dropna()
        X = combo[[c]]
        y = combo[['pheno']]
        
        logit_model=sm.Logit(y,X)
        result=logit_model.fit()
        auc = roc_auc_score(y, result.predict())
        pval = result.pvalues[0]
        
        
        for i in range(0,3):
            X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.1, random_state=random_states[i])
            logit_model=sm.Logit(y_train,X_train)
            result=logit_model.fit()
            if len(y_test['pheno'].unique())>1: # check to be sure we have more than 1 pheno in y_test
                auc_split = roc_auc_score(y_test, result.predict(X_test))
                aurocs.append(auc_split)
            
        aurocs = pd.Series(aurocs)
        
        logreg_out[c]={'LogReg P-value':pval, 'AUROC (Logistic Regression)':auc, 'LogReg AUROC Lower Limit (95%CI)':aurocs.min(), 'LogReg AUROC Upper Limit (95%CI)':aurocs.max()}
        
    logreg_out = pd.DataFrame(logreg_out).transpose()
    print('TODO add CI/pval?? - see Bob 2/28/21 email')
        # embed()
    return logreg_out
    
    ''' END NEW LOGREG SECTION '''
    


def get_outliers(controlsC, samplesC):
    def reject_outliers(data, m = 2.):
        d = np.abs(data - np.median(data))
        # embed()
        mdev = np.median(d)
        s = d/mdev if mdev else 0.
        min = np.median(data)-m*mdev
        max = m*mdev+np.median(data)
        try:
            return data[s<m], min, max
        except:
            print('error on outlier rejection')
            return data, min, max
    
    def get_outliers(data,name):
        # this comes from the post here: https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
        outliers = {}
        for c in data.columns:
            res, min, max = reject_outliers(data[c].dropna(),m=4)
            # print c, len(data)-len(res)
            num_outliers = len(data)-len(res)
            if num_outliers==20:
                num_outliers=0
            if num_outliers>0:
                outliers[c] = {'%s Outliers'%name:'YES','%s num_outliers'%name:num_outliers,'%s Low Outlier Threshold'%name:min,'%s High Outlier Threshold'%name:max}
            else:
                outliers[c] = {'%s Outliers'%name:'NO','%s num_outliers'%name:num_outliers,'%s Low Outlier Threshold'%name:min,'%s High Outlier Threshold'%name:max}
        
        return pd.DataFrame(outliers).transpose()
    
    control_outliers = get_outliers(controlsC,'Control')
    sample_outliers = get_outliers(samplesC,'Case')
    outliers = pd.merge(control_outliers, sample_outliers, left_index=True, right_index=True).dropna()
    
    return outliers  

def calc_other_stats(controls, controls_linear, controlsC, samples, samples_linear, samplesC):
    
    # Mann Whitney https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
    from scipy import stats 
    
    mann_whitney = {}
    mann_whitney2 = {}
    # embed()
    for c in controlsC.columns:
        stat, p_mannw2 = stats.mannwhitneyu(controlsC[c].dropna(),samplesC[c].dropna(), use_continuity=False, alternative='two-sided') # NOTE: these parameters will make the MWpval=KWpval (for two samples)
        stat, p_mannw = stats.mannwhitneyu(controlsC[c].dropna(),samplesC[c].dropna()) # NOTE: these parameters will make the MWpval=KWpval (for two samples)
        mann_whitney[c]=p_mannw
        mann_whitney2[c]=p_mannw2
    mann_whitney = pd.Series(mann_whitney,name='Mann-Whitney p (1-sided)')
    mann_whitney2 = pd.Series(mann_whitney2,name='Mann-Whitney U test p')
    
    # normal distr:
    # D'agostino https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html
    dagostino = {}
    # embed()
    for c in controlsC.columns:
        if len(controlsC[c].dropna())==0:
            continue
        if len(controlsC[c].dropna())>8:
            k2, p_dagostino = stats.normaltest(controlsC[c].dropna())
            dagostino[c] = p_dagostino
        else:
            dagostino[c] = 'Too few samples'
    dagostino = pd.Series(dagostino, name='D\'Agostino & Pearson omnibus p')
    
    dagostino_linear = {}
    
    for c in controls.columns[2:]:
        if len(controlsC[c].dropna())==0:
            continue
        elif len(controlsC[c].dropna())>8:
            k2, p_dagostino = stats.normaltest(pd.to_numeric(controls[c]).dropna())
            dagostino_linear[c] = p_dagostino
        else:
            dagostino_linear[c] = 'Too few samples'
    dagostino_linear = pd.Series(dagostino_linear, name='D\'Agostino & Pearson omnibus p (Linear)')
    
    
    # Shapiro & Wilk https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html
    shapiro_wilk = {}
    for c in controlsC.columns:
        if len(controlsC[c].dropna())==0:
            continue
        w, p_shapiro = stats.shapiro(controlsC[c].dropna())
        shapiro_wilk[c]=p_shapiro
    shapiro_wilk = pd.Series(shapiro_wilk, name='Shapiro-Wilk normality p')
    
    shapiro_wilk_linear = {}
    for c in controls.columns[2:]:
        if len(controls[c].dropna())==0:
            continue
        w, p_shapiro = stats.shapiro(controls[c].dropna())
        shapiro_wilk_linear[c]=p_shapiro
    shapiro_wilk_linear = pd.Series(shapiro_wilk_linear, name='Shapiro-Wilk normality p (Linear)')
    
    out = pd.DataFrame([mann_whitney,mann_whitney2,dagostino,shapiro_wilk,shapiro_wilk_linear,dagostino_linear]).transpose().dropna()
    
    return out
    

    
def calc_sem_etc(controls_linear, controlsC, samples_linear, samplesC):
    from scipy import stats
    ''' ADDITIONS TO TABLE S1 (SEM, SD 1.96 ETC) 6/12/2018 '''
    # followed work in D:\Dropbox (Personal)\naviauxlab_informatics\Dutch_Shepherd\analyses_for_updated_tableS1\Dog_No154141_final_with replicates_noLipoic (1).xlsx
    
    #log space:
    control_log_mean = controlsC.mean()
    control_log_mean.name = 'Control Mean'
    
    control_log_std = controlsC.std()
    control_log_std_star = 2**controlsC.std()
    control_log_std_star.name = 'Control SD*'
    
    
    control_sem = controlsC.sem()
    control_sem_star = 2**controlsC.sem()
    control_sem_star.name = 'Control SEM*'
    
    case_log_std = samplesC.std()
    case_log_std_star = 2**samplesC.std()
    case_log_std_star.name = 'Case SD*'

    case_sem = samplesC.sem()
    case_sem_star = 2**samplesC.sem()
    case_sem_star.name = 'Case SEM*'
    
    case_log_mean = samplesC.mean()
    case_log_mean.name = 'Case Mean'
    
    
    # calculating standard error of means 9/19/2018
    sda = controlsC.std()**2
    sdb = samplesC.std()**2
    
    sem = np.sqrt(sda/len(controlsC)+sdb/len(samplesC))
    sem.name = 'Log2 SEM'
    lsem = 2**sem
    lsem.name = 'Linear SEM'
    
    tmp = controlsC.append(samplesC)
    
    
    # embed()
    

    control_linear_mean = 2**control_log_mean # geoMean
    control_linear_mean.name='Linear Control Mean (Geomean)'
    
    linear_std = 2**control_log_std
    linear_std.name='Linear SD*'
    
    linear_sem = 2**control_sem
    linear_sem.name = 'Linear SEM*'
   
    
    # control_log_std_196 = control_log_std**1.96
    # control_log_std_196.name = 'SD*'
    control_low25 = control_linear_mean/(control_log_std_star**1.96)
    control_low25.name = '2.50%'
    control_high25 = control_linear_mean*(control_log_std_star**1.96)
    control_high25.name = '97.50%'
    
    # HiLoRatio
    hilo = control_high25/control_low25
    hilo.name = 'Hi/Lo Ratio'
    
    #Linear Case Mean (Geomean)
    case_linear_mean = 2**case_log_mean
    case_linear_mean.name = 'Linear Case Mean (Geomean)'
    # embed()
    
    
    LogLinFC = case_linear_mean/control_linear_mean
    LogLinFC.name = 'Fold Change'
    
    pvals2 = {}
    for c in controls_linear.columns:
        res, pval = stats.ttest_ind(controls_linear[c].dropna(), samples_linear[c].dropna())
        pvals2[c]=pval
    calc_pvals = pd.Series(pvals2)
    calc_pvals.name = 'P Value - Linear'
    
    controlsMedian = controlsC.median()
    controlsMedian.name = 'Control Median'
    samplesMedian = samplesC.median()
    samplesMedian.name = 'Case Median'
    
    out = pd.DataFrame([sem, lsem, control_log_mean,control_log_std_star,control_sem_star,case_log_std_star, case_sem_star,control_low25,control_high25,case_log_mean,LogLinFC,calc_pvals,hilo,linear_sem,linear_std,control_linear_mean,controlsMedian,case_linear_mean,samplesMedian]).transpose()
    
    return out.dropna()
    
    
def calc_stats(samples, samples_linear, samplesC, controls, controls_linear, controlsC):
    from scipy import stats
    overlap = list(set(controlsC.columns&samplesC.columns))
    overlap.sort()

    # linear z score
    
    linear_zscore = pd.DataFrame(samples_linear.mean().sub(controls_linear.mean()).div(controls_linear.std()).reset_index().values.copy(), columns=['MRM Name', 'Linear Z Score']).dropna()
    linear_zscore.index = linear_zscore['MRM Name']
    
    # calc t tests:
    
    
    len_control = len(controlsC)
    len_samples = len(samplesC)
    filter_samples = []
    out2 = {}
    out3 = {}
    out4 = {}
    for o in overlap:
        controlC = controlsC[o].dropna()
        sampleC = samplesC[o].dropna()
        if len(controlC)>len_control*0.7 and len(sampleC)>len_samples*0.7:
        
            rej, pval = stats.ttest_ind(controlC, sampleC)
            out2[o]=pval
            stat, kw_pval = stats.kruskal(controlC, sampleC)
            # stat2, kw_pval2 = scipy.stats.mstats.kruskalwallis(controlsC[o], samplesC[o])
            
            
            out3[o]=kw_pval
            rej, welch_pval = stats.ttest_ind(controlC, sampleC, equal_var=False)
            out4[o]=welch_pval
            
        else:
            # print 'filter',o,'too few datapoints'
            filter_samples.append(o)
    
    pvals = pd.DataFrame(pd.Series(out2).dropna(),columns=['P Value'])
    kw_pvals = pd.DataFrame(pd.Series(out3).dropna(),columns=['Kruskal-Wallis p'])
    welch_pvals = pd.DataFrame(pd.Series(out4).dropna(),columns=['Welch\'s ttest p'])
    welch_pvals = pd.DataFrame(pd.Series(out4).dropna(),columns=['Welch\'s ttest p'])
    
    FC = pd.DataFrame(samples.mean().dropna()/controls.mean().dropna(),columns=['Fold Change_log'])
    if len(FC)==0:
        FC = pd.DataFrame(samples[samples.columns[2:]].mean().dropna()/controls[controls.columns[2:]].mean().dropna(),columns=['Fold Change_log'])
    
    out = pd.merge(FC, pvals, left_index=True, right_index=True)
    out = pd.merge(out, kw_pvals, left_index=True, right_index=True)
    out = pd.merge(out, welch_pvals, left_index=True, right_index=True)
    out = pd.merge(out, linear_zscore, left_index=True, right_index=True)
    
    del out['MRM Name']
    
    num_mets_measured = len(out)
    
    out.sort_values('P Value', inplace=True)
    # out = out.reset_index()
    # out = pd.merge(out, db, left_on='index', right_on='MRM Name_y', how='left')
    

    #from http://jpktd.blogspot.com/2013/04/multiple-testing-p-value-corrections-in.html
    import statsmodels.stats.multitest as smm
    from statsmodels.stats.multitest import fdrcorrection
    from statsmodels.stats.multitest import fdrcorrection_twostage

    def FDR(p):
        
        """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
        p = np.asfarray(p)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / np.arange(len(p), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]
    
    pvals=out['P Value']
    
    res = fdrcorrection_twostage(pvals, alpha=0.05, method='bh')
    rej = res[0]
    pval_corr = res[1]
    
    fdr = pd.DataFrame(pd.Series(pval_corr), columns=['FDR'])
    fdr.index=pvals.index
    out = pd.merge(out, fdr, left_index=True, right_index=True)
    
    ''' calculating storey's qvalue  '''
    # uses package: https://github.com/nfusi/qvalue'''
    from . import qvalue
    # embed()
    qv = qvalue.estimate(pvals) # Storey's q value
    qv = pd.DataFrame(pd.Series(qv), columns=['q-Value'])
    qv.index=pvals.index
    out = pd.merge(out, qv, left_index=True, right_index=True)
    
    return out
    


def calc_heteroscasdicity(controls_linear, controlsC, samples_linear, samplesC):
    ''' ADDING HET_BREUSCHPAGAN TEST FOR heteroscedasticity  12/22/2018 '''
    # https://www.statsmodels.org/dev/generated/statsmodels.stats.diagnostic.het_breuschpagan.html
    het_b_out={}
    import statsmodels.api as sm
    from statsmodels.stats.diagnostic import het_breuschpagan
    
    # 1. Control Heteroscedasticity (log space Breusch-Pagan p)
    # 2. Control Heteroscedasticity (linear space Breusch-Pagan p)
    # 3. Case Heteroscedasticity (log space Breusch-Pagan p)
    # 4. Case Heteroscedasticity (linear space Breusch-Pagan p)
    data_sources = {
        'Control Heteroscedasticity (linear space Breusch-Pagan p)':controls_linear,
        'Control Heteroscedasticity (log space Breusch-Pagan p)':controlsC,
        'Case Heteroscedasticity (linear space Breusch-Pagan p)':samples_linear,
        'Case Heteroscedasticity (log space Breusch-Pagan p)':samplesC
    }
    # embed()
    

    print('running het_breuschpagan')
    
    
    for c in controlsC:
        
        # tesing if combining control/case changes results
        exog = data[c].dropna()
        if len(exog)>0:
            endog = data['Group']
            endog[endog==exp]=1
            endog[endog==control]=0
            sm.add_constant(exog, prepend=False) # not sure whether to add this or not...
            # embed()
            mod = sm.OLS(pd.to_numeric(endog), exog)
            res = mod.fit()
            lm, lm_pvalue,fvalue, f_pvalue = het_breuschpagan(res.resid,np.array([exog]).transpose())
            # embed()
            het_b_out[c]={}
            for source in data_sources.keys():
                data_source = data_sources[source]
                # c = controls.columns[2]
                # tmp = sm.add_constant(controls[c], prepend=False)
                endog = [1]*len(data_source)
                exog = data_source[c]
                mod = sm.OLS(endog, exog)
                
                res = mod.fit()
                
                # print(res.summary())
                lm, lm_pvalue,fvalue, f_pvalue = het_breuschpagan(res.resid,np.array([exog]).transpose())
                # print c, f_pvalue
                het_b_out[c][source]=f_pvalue
    het_b_out = pd.DataFrame(het_b_out).transpose()



def build_table_s1(outfile, study_name, data, control, exp, table1_mets, db_file, ML_iter=10, VIP_cutoff='', mda_file='', run_het=False,       knn_vip_threshold=1):
    out = {}
    data_file = data
    data = pd.read_csv(data)
    
    num_mets_targeted = len(data.columns)-2 # assuming first two are "Sample" and "Group"
    
    mrm = pd.read_excel(db_file,'MRM')
    alias = pd.read_excel(db_file,'Alias')
    chem = pd.read_excel(db_file,'Chemical')
    db = pd.merge(mrm,alias, left_on='MRM Index', right_on='MRM Index')[['MRM Name_x','MRM Index','MRM Name_y','Chemical Index']]
    db = pd.merge(db,chem[['Chemical Index','CAS Number','HMDB number']], left_on='Chemical Index', right_on='Chemical Index', how='left')
    # embed()
    
    controls=data[data.Group==control]
    samples=data[data.Group==exp]

    # calc t tests:
    from scipy import stats
    controls_linear = controls[controls.columns[2:]]
    controlsC = controls[controls.columns[2:]].apply(np.log2).replace([np.inf,-np.inf],np.nan)
    # embed()
    samplesC = samples[samples.columns[2:]].apply(np.log2).replace([np.inf,-np.inf],np.nan)
    samples_linear = samples[samples.columns[2:]]
    overlap = list(set(controlsC.columns&samplesC.columns))
    overlap.sort()
    
    # linear z score
    
    linear_zscore = pd.DataFrame(samples_linear.mean().sub(controls_linear.mean()).div(controls_linear.std()).reset_index().values.copy(), columns=['MRM Name', 'Linear Z Score']).dropna()
    linear_zscore.index = linear_zscore['MRM Name']
    # embed()
    
    len_control = len(controlsC)
    len_samples = len(samplesC)
    filter_samples = []
    # embed()
    out2 = {}
    out3 = {}
    out4 = {}
    for o in overlap:
        controlC = controlsC[o].dropna()
        sampleC = samplesC[o].dropna()
        if len(controlC)>len_control*0.7 and len(sampleC)>len_samples*0.7:
        
            rej, pval = stats.ttest_ind(controlC, sampleC)
            out2[o]=pval
            stat, kw_pval = stats.kruskal(controlC, sampleC)
            # stat2, kw_pval2 = scipy.stats.mstats.kruskalwallis(controlsC[o], samplesC[o])
            
            
            out3[o]=kw_pval
            rej, welch_pval = stats.ttest_ind(controlC, sampleC, equal_var=False)
            out4[o]=welch_pval
            
        else:
            print('filter',o,'too few datapoints')
            filter_samples.append(o)
        # embed()
    pvals = pd.DataFrame(pd.Series(out2).dropna(),columns=['P Value'])
    kw_pvals = pd.DataFrame(pd.Series(out3).dropna(),columns=['Kruskal-Wallis p'])
    welch_pvals = pd.DataFrame(pd.Series(out4).dropna(),columns=['Welch\'s ttest p'])
    welch_pvals = pd.DataFrame(pd.Series(out4).dropna(),columns=['Welch\'s ttest p'])
    # embed()
    # test1 = controlsC['SM(d18:1/26:0 OH)']
    # test2 = samplesC['SM(d18:1/26:0 OH)']
    # test = stats.ttest_ind(test1, test2, equal_var=True)
    # embed()
    # embed() # TO DO determine why excel ttest gives different results from scipy...
    # pvals = stats.ttest_ind(controlsC[overlap], samplesC[overlap])
    
    # embed()
    
    # pvals = pd.DataFrame(pd.Series(test.pvalue, index=overlap),columns=['pval']) # new version of numpy/stats
    # pvals = pd.DataFrame(pd.Series(pvals[1], index=overlap),columns=['P Value']).dropna()
    # pvals_copy = pvals
    # pvals = pd.merge(pvals, maleC.transpose(), left_index=True, right_index=True, how='left')
    # pvals = pd.merge(pvals, femaleC.transpose(), left_index=True, right_index=True, how='left', suffixes=['_male','_female'])
    # pvals.to_csv('test.csv')

    FC = pd.DataFrame(samples.mean().dropna()/controls.mean().dropna(),columns=['Fold Change_log'])
    # linear_zscore = linear_zscore.dropna()
    # embed()
    out = pd.merge(FC, pvals, left_index=True, right_index=True)
    out = pd.merge(out, kw_pvals, left_index=True, right_index=True)
    out = pd.merge(out, welch_pvals, left_index=True, right_index=True)
    out = pd.merge(out, linear_zscore, left_index=True, right_index=True)
    # embed()
    del out['MRM Name']
    # embed()
    num_mets_measured = len(out)
    # out = pd.merge(FC, pvals)
    out.sort_values('P Value', inplace=True)
    out = out.reset_index()
    # embed()out
    # embed()
    # out = pd.merge(out, db, left_on='MRM Name', right_on='MRM Name_y', how='left')
    out = pd.merge(out, db, left_on='index', right_on='MRM Name_y', how='left')
    # out = out.rename(columns={'MRM Name_x':'Study MRM Name'})
    # out = pd.merge(out, db, left_on='index', right_on='MRM Name_y', how='left') # edited 1/13/2021 to deal with new PTLDS
    
    
    # check to be sure all mets map over
    tmp = out.fillna('-')
    tmp = tmp[tmp['MRM Name_x']=='-'] #Study MRM Name
    if len(tmp)>0:
        print('dB mapping error')
        print(tmp)
        # embed()
    
    # embed()

    #from http://jpktd.blogspot.com/2013/04/multiple-testing-p-value-corrections-in.html
    import statsmodels.stats.multitest as smm
    from statsmodels.stats.multitest import fdrcorrection
    from statsmodels.stats.multitest import fdrcorrection_twostage

    def FDR(p):
        
        """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
        p = np.asfarray(p)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / np.arange(len(p), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]
    
    pvals=out['P Value']
    # embed()
    
    ''' calculating adjusted p-vals, e.g FDR '''
    
    # test = FDR(pvals.values) # this gives same value as below:
    
    # alpha = 0.1
    # method = 'fdr_i'
    # rej, pval_corr = smm.multipletests(pval_raw, alpha=alpha, method=method)
    # rej, pval_corr = smm.multipletests(pvals, method='fdr_bh')[:2]
    # rej, pvalscorr, ntests, alpha_stages = fdrcorrection_twostage(pvals, alpha=alpha)
    
    # rej, pval_corr = fdrcorrection(pvals, alpha=0.05)
    res = fdrcorrection_twostage(pvals, alpha=0.05, method='bh')
    rej = res[0]
    pval_corr = res[1]
    
    # embed()
    
    # import pylab as plt
    # pvals.hist()
    # plt.show()
    
    fdr = pd.DataFrame(pd.Series(pval_corr), columns=['FDR'])
    
    out = pd.merge(out, fdr, left_index=True, right_index=True)
    
    # print 'WHY are pvals so HIGH?' # TODO
    # embed()
    
    
    
    ''' calculating storey's qvalue  '''
    # uses package: https://github.com/nfusi/qvalue'''
    from . import qvalue
    qv = qvalue.estimate(pvals) # Storey's q value
    
    qv = pd.DataFrame(pd.Series(qv), columns=['q-Value'])
    # embed()
    out = pd.merge(out, qv, left_index=True, right_index=True)
    
    # BELOW CALCULATING POWER REMOVED 12/12/2017 for CFS2 results
    ''' calculating POWER '''
    '''
    # getting this from here: http://jpktd.blogspot.com/2013/03/statistical-power-in-statsmodels.html
    import statsmodels.stats.power as smp
    # not sure this is correct...
    # power = 1-smp.TTestIndPower().solve_power(pvals.values, nobs1=30, ratio=1, alpha=0.05, alternative='two-sided')
    # test2 = smp.ttest_power(pvals, nobs=60, alpha=0.1, alternative='two-sided')
    # ABOVE doesnt work
    
    # TEST calcing my own power
    def calc_power(control,sample):
        import math
        from scipy import stats
        s1 = control
        s2 = sample
        
        # if control.name=='FAD':
            # import pylab as plt
            # s1.hist(bins=30)
            # s2.hist(bins=30)
            # plt.show()
            
        # embed()
        sem = s2.std()/math.sqrt(len(s2))
        crit_val = stats.norm.ppf(0.99, loc=s2.mean(), scale=sem) # JM changed to 0.99 from 0.95

        s1_sem = s1.std()/math.sqrt(len(s1))
        beta = stats.norm.cdf(crit_val, s1.mean(), s1_sem)
        # if control.name=='FAD':
            # embed()
        # print met, 1-beta
        power = 1-beta
        # embed()
        return power
        
    pcontrol = controlsC[overlap]
    psample = samplesC[overlap]
    out['Power'] = ''
    for i in out.index:
        met = out.loc[i,'index']
        s1 = pcontrol[met]
        s2 = psample[met]
        # print met
        # print ttest_ind(s1,s2)
        # power = calc_power(s1,s2)
        if s1.mean()>s2.mean():
            power = calc_power(s1,s2)
        else:
            power = calc_power(s2,s1)
        out.loc[i,'Power']=power
    # embed()
    # adding by index above so no need to do merge below:
    # power = pd.DataFrame(pd.Series(power), columns=['Power'])
    # out = pd.merge(out, power, left_index=True, right_index=True)
    '''
    
    # BELOW RANDOM FOREST AND MDA REMOVED ON 12/12/2017 for CFS2 figs
    ''' CALCULATING Random forest out of box for each met '''
    ''' calculate RFOOB for each met '''
    '''
    from sklearn.ensemble import RandomForestClassifier
    
    from sklearn.model_selection import train_test_split
    # from sklearn.datasets import make_classification
    # X, y = make_classification(n_samples=1000, n_features=4,
                           # n_informative=2, n_redundant=0,
                           # random_state=0, shuffle=False)
    clf = RandomForestClassifier(max_depth=2, random_state=0)
    
    rf_control = controlsC[overlap]
    rf_sample = samplesC[overlap]
    
    X_full = pd.concat([rf_sample,rf_control])
    X_full = X_full.transpose().dropna().transpose()
    Y = pd.Series([1]*len(rf_sample)+[0]*len(rf_control))
    
    rfoob = {}
   
    
    def get_best_features(X,Y, num_features=7):
        sample = X.ix[Y[Y==1].index]
        control = X.ix[Y[Y==0].index] 
        res = stats.ttest_ind(control, sample)
        out = {}
        for c in sample.columns:
            rej, pval = stats.ttest_ind(control[c],sample[c])
            # print c, pval
            out[c] = pval
        out = pd.Series(out)
        out.sort()
        return out
    bpv = get_best_features(X_full, Y)
    
    out = pd.merge(out,pd.DataFrame(bpv,columns=['BPV']),left_on='MRM Name_y', right_index=True)
    # embed()
    
    # embed()
    
    feature_count = 7
    top_features = bpv[:feature_count].index
    # RUN ONCE ON ALL
    score_list = []
    for ii in range(0,ML_iter):
        X = X_full[top_features]
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=i)
        train_matrix = X_train
        test_matrix = X_test
        train_label = y_train
        true_test_label = y_test
        
               
        regr =RandomForestClassifier()
        model_fit=regr.fit(train_matrix,train_label)
        model_score=regr.score(test_matrix,true_test_label)
        predict_label=regr.predict(test_matrix)
        
        # print model_score
        score_list.append(model_score)
    avg_score = sum(score_list)/len(score_list)
    mean_acc = avg_score
    print 'MEAN accuracy with all mets is', mean_acc
    
    
    # updated to REMOVE each met at a time
    for i in range(0,len(X_full.columns)):
        score_list = []
        c = X_full.columns[i]
        print 'Running ML for',c , ML_iter, 'iterations' 
        for ii in range(0,ML_iter):
            
            if c in top_features:
                tmp_top_features = bpv[:feature_count+1]
                tmp = tmp_top_features.index
                tmp_top_features = tmp[tmp!=c]
                X = X_full[tmp_top_features]
                # embed()
            else:
                
                X = X_full[top_features]
            # embed()    
            # X = X_full[X_full.columns[i:i+1]]
            # X1 = X_full[X_full.columns[:i]]
            # X2 = X_full[X_full.columns[i+1:]]
            # X = pd.concat([X1,X2],axis=1) # REMOVE the met in question
            
            
            X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=i)
            train_matrix = X_train
            test_matrix = X_test
            train_label = y_train
            true_test_label = y_test
            
                   
            regr =RandomForestClassifier()
            model_fit=regr.fit(train_matrix,train_label)
            model_score=regr.score(test_matrix,true_test_label)
            predict_label=regr.predict(test_matrix)
            
            # print model_score
            score_list.append(model_score)
        avg_score = sum(score_list)/len(score_list)
        MDA = abs(mean_acc-avg_score)
        print '\tMDA = ', MDA
        # print X.columns
        # if c in top_features:
            # embed()
        rfoob[i] = MDA
    rfoob = pd.DataFrame(pd.Series(rfoob), columns=['MDA'])
    # embed()
    out = pd.merge(out, rfoob, left_index=True, right_index=True)
    '''
    
    
    
    # out = pd.merge(out,het_b_out, left_on='MRM Name_y', right_index=True)    
    # embed()
       
    
    ''' ADDITIONS TO TABLE S1 (SEM, SD 1.96 ETC) 6/12/2018 '''
    # followed work in D:\Dropbox (Personal)\naviauxlab_informatics\Dutch_Shepherd\analyses_for_updated_tableS1\Dog_No154141_final_with replicates_noLipoic (1).xlsx
    # embed()
    # controls
    # for testing: 
    met = '1-Methyladenosine'
    
    #log space:
    control_log_mean = controlsC.mean()
    control_log_mean.name = 'Control Mean'
    
    control_log_std = controlsC.std()
    control_log_std_star = 2**controlsC.std()
    control_log_std_star.name = 'Control SD*'
    
    
    control_sem = controlsC.sem()
    control_sem_star = 2**controlsC.sem()
    control_sem_star.name = 'Control SEM*'
    
    case_log_std = samplesC.std()
    case_log_std_star = 2**samplesC.std()
    case_log_std_star.name = 'Case SD*'

    case_sem = samplesC.sem()
    case_sem_star = 2**samplesC.sem()
    case_sem_star.name = 'Case SEM*'
    
    case_log_mean = samplesC.mean()
    case_log_mean.name = 'Case Mean'
    
    
    # calculating standard error of means 9/19/2018
    sda = controlsC.std()**2
    sdb = samplesC.std()**2
    
    sem = np.sqrt(sda/len(controlsC)+sdb/len(samplesC))
    sem.name = 'Log2 SEM'
    lsem = 2**sem
    lsem.name = 'Linear SEM'
    
    tmp = controlsC.append(samplesC)
    
    
    # embed()
    

    control_linear_mean = 2**control_log_mean # geoMean
    control_linear_mean.name='Linear Control Mean (Geomean)'
    
    linear_std = 2**control_log_std
    linear_std.name='Linear SD*'
    
    linear_sem = 2**control_sem
    linear_sem.name = 'Linear SEM*'
   
    
    # control_log_std_196 = control_log_std**1.96
    # control_log_std_196.name = 'SD*'
    control_low25 = control_linear_mean/(control_log_std_star**1.96)
    control_low25.name = '2.50%'
    control_high25 = control_linear_mean*(control_log_std_star**1.96)
    control_high25.name = '97.50%'
    
    # HiLoRatio
    hilo = control_high25/control_low25
    hilo.name = 'Hi/Lo Ratio'
    
    #Linear Case Mean (Geomean)
    case_linear_mean = 2**case_log_mean
    case_linear_mean.name = 'Linear Case Mean (Geomean)'
    # embed()
    
    
    LogLinFC = case_linear_mean/control_linear_mean
    LogLinFC.name = 'Fold Change'
    # rejs2, pvals2 = stats.ttest_ind(controls, samples) # edit for two-tailed, two-sample equal variants, homeoscedastic 
    rejs2, pvals2 = stats.ttest_ind(controls_linear, samples_linear) # edit for two-tailed, two-sample equal variants, homeoscedastic 
    
    
    calc_pvals = pd.Series(pvals2,index=controlsC.columns)
    calc_pvals.name = 'P Value - Linear'
    
    # Mann Whitney https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
    from scipy import stats 
    
    mann_whitney = {}
    mann_whitney2 = {}
    # embed()
    for c in controlsC.columns:
        stat, p_mannw2 = stats.mannwhitneyu(controlsC[c].dropna(),samplesC[c].dropna(), use_continuity=False, alternative='two-sided') # NOTE: these parameters will make the MWpval=KWpval (for two samples)
        stat, p_mannw = stats.mannwhitneyu(controlsC[c].dropna(),samplesC[c].dropna()) # NOTE: these parameters will make the MWpval=KWpval (for two samples)
        mann_whitney[c]=p_mannw
        mann_whitney2[c]=p_mannw2
    mann_whitney = pd.Series(mann_whitney,name='Mann-Whitney p (1-sided)')
    mann_whitney2 = pd.Series(mann_whitney2,name='Mann-Whitney U test p')
    
        
    # normal distr:
    # D'agostino https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html
    dagostino = {}
    # embed()
    for c in controlsC.columns:
        if len(controlsC[c])>8:
            k2, p_dagostino = stats.normaltest(controlsC[c])
            dagostino[c] = p_dagostino
        else:
            dagostino[c] = 'Too few samples'
    dagostino = pd.Series(dagostino, name='D\'Agostino & Pearson omnibus p')
    
    # embed()
    
    dagostino_linear = {}
    for c in controls.columns[2:]:
        if len(controlsC[c])>8:
            k2, p_dagostino = stats.normaltest(controls[c])
            dagostino_linear[c] = p_dagostino
        else:
            dagostino_linear[c] = 'Too few samples'
    dagostino_linear = pd.Series(dagostino_linear, name='D\'Agostino & Pearson omnibus p (Linear)')
    
    
    # Shapiro & Wilk https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html
    shapiro_wilk = {}
    for c in controlsC.columns:
        w, p_shapiro = stats.shapiro(controlsC[c])
        shapiro_wilk[c]=p_shapiro
    shapiro_wilk = pd.Series(shapiro_wilk, name='Shapiro-Wilk normality p')
    
    shapiro_wilk_linear = {}
    for c in controls.columns[2:]:
        w, p_shapiro = stats.shapiro(controls[c])
        shapiro_wilk_linear[c]=p_shapiro
    shapiro_wilk_linear = pd.Series(shapiro_wilk_linear, name='Shapiro-Wilk normality p (Linear)')
    
    def reject_outliers(data, m = 2.):
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d/mdev if mdev else 0.
        min = np.median(data)-m*mdev
        max = m*mdev+np.median(data)
        
        return data[s<m], min, max
    # trying with pierce: trouble...
    # n=1
    # m=1
    # res = peirce_dev(20, n, m)
    
    ''' ADDING log reg 2/27/21 '''
    
    tmp_controlsC = controlsC.copy()
    tmp_controlsC['pheno'] = 0
    tmp_samplesC = samplesC.copy()
    tmp_samplesC['pheno'] = 1

    logreg_data = tmp_controlsC.append(tmp_samplesC)
    logreg_cnts = logreg_data.count()
    to_keep = logreg_cnts[logreg_cnts>0]
    logreg_data = logreg_data[to_keep.index]
    
    import statsmodels.api as sm
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score
    logreg_out = []
    
    for c in logreg_data.columns[:-1]:
        X = logreg_data[[c]]
        y = logreg_data[['pheno']]
        
        X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.33, random_state=42)
        
        logit_model=sm.Logit(y,X)
        result=logit_model.fit()
        auc = roc_auc_score(y, result.predict())
        pval = result.pvalues[0]
        logit_model=sm.Logit(y_train,X_train)
        result=logit_model.fit()
        auc_split = roc_auc_score(y_test, result.predict(X_test))
        pval_split = result.pvalues[0]
        logreg_out.append({'LR MRM Name':c, 'LogReg_p':pval, 'LogReg_AUC':auc, 'LogReg_p(split)':pval_split, 'LogReg_AUC(split)':auc_split})
    logreg_out = pd.DataFrame(logreg_out)
    # embed()
    
    ''' END NEW LOGREG SECTION '''
    
    
    def get_outliers(data,name):
        # this comes from the post here: https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
        outliers = {}
        for c in data.columns:
            res, min, max = reject_outliers(data[c],m=4)
            # print c, len(data)-len(res)
            num_outliers = len(data)-len(res)
            if num_outliers==20:
                num_outliers=0
            if num_outliers>0:
                outliers[c] = {'%s Outliers'%name:'YES','num_outliers':num_outliers,'%s Low Outlier Threshold'%name:min,'%s High Outlier Threshold'%name:max}
            else:
                outliers[c] = {'%s Outliers'%name:'NO','num_outliers':num_outliers,'%s Low Outlier Threshold'%name:min,'%s High Outlier Threshold'%name:max}
        
        return pd.DataFrame(outliers).transpose()
    
    control_outliers = get_outliers(controlsC,'Control')
    sample_outliers = get_outliers(samplesC,'Case')
    
    controlsMedian = controlsC.median()
    controlsMedian.name = 'Control Median'
    samplesMedian = samplesC.median()
    samplesMedian.name = 'Case Median'
    # embed()
    
    # updated 9/19/2018
    new_analysis = pd.DataFrame([sem, lsem, control_log_mean,control_log_std_star,control_sem_star,case_log_std_star, case_sem_star,control_low25,control_high25,case_log_mean,LogLinFC,calc_pvals,hilo,linear_sem,linear_std,control_linear_mean,controlsMedian,case_linear_mean,samplesMedian,mann_whitney,mann_whitney2,dagostino,shapiro_wilk,shapiro_wilk_linear,dagostino_linear]).transpose()
    # new_analysis = pd.DataFrame([control_log_mean,control_log_std_star,control_sem_star,case_log_std_star, case_sem_star,control_low25,control_high25,case_log_mean,LogLinFC,calc_pvals,hilo,linear_sem,linear_std,control_linear_mean,controlsMedian,case_linear_mean,samplesMedian,mann_whitney,dagostino,shapiro_wilk,shapiro_wilk_linear,dagostino_linear]).transpose()
    # embed()
    new_analysis = pd.merge(new_analysis, control_outliers, left_index=True, right_index=True)
    new_analysis = pd.merge(new_analysis, sample_outliers, left_index=True, right_index=True)
    # embed()
    
    # metaboanalysts changes met names, use below to fix
    new_analysis['metaboanalyst_name']=new_analysis.index
    new_analysis['metaboanalyst_name'] = new_analysis['metaboanalyst_name'].str.replace(':','')
    new_analysis['metaboanalyst_name'] = new_analysis['metaboanalyst_name'].str.replace('(','')
    new_analysis['metaboanalyst_name'] = new_analysis['metaboanalyst_name'].str.replace(')','')
    
    ''' MEAN DECREASE IN ACCURACY '''
    # TOO SLOW NOW in PYTHON
    # USE results from metaboanalyst for now
    if mda_file!='':
        mda = pd.read_csv(mda_file, index_col=0)
        mda = mda.rename(columns={'MeanDecreaseAccuracy':'RFA MDA (5000 trees)'})
        
        new_analysis = pd.merge(new_analysis, mda, left_on='metaboanalyst_name', right_index=True, how='left')
    
    ''' KMEANS ADDITIONS '''
        
    from sklearn.cluster import KMeans
    
    
    
    data_col=2
    
    data2 = {}
    for c in data.columns[data_col:]:
        try:
            data2[c] = pd.to_numeric(data[c])
        except:
            print('error', c)
            embed()
    data2 = pd.DataFrame(data2, index=data.index)
    data2 = np.log(data2)
    
    
    orig_data = pd.read_csv(data_file, index_col=0)
    
    # updated 10/6/2020 for zscore tab
    zscores = get_zscores(orig_data, exp, control)
    zscores = zscores['case'].append(zscores['control'])
    zscores['Group'] = orig_data['Group']
    zscores = zscores[['Group'] + zscores.columns[:-1].tolist()]
    
    #updated data3 to come from these zscore rather than using scipy:
    
    # OLD CODE:
    # data3 = stats.zscore(data2,axis=0) # axis=0 is columnwise
    
    #NEW CODE:
    data3 = zscores
    
    data3 = pd.DataFrame(data3,columns=data2.columns)
    test=controlsC.append(samplesC) # THIS IS ASSUMING  use kmeans on ALL DATA
    # test=controlsC # THIS IS ASSUMING  use kmeans on controls
    # test=samplesC # THIS IS ASSUMING  use kmeans on samples
    
    # 10/5/2020: errors in knn if any analytes have N/A - e.g. FGF in Males has errors
    # cnts = data3.count()
    # to_keep = cnts[cnts>0]
    
    X = data3.transpose().dropna()
    
    
    # embed()
    kmeans_clusters = [5, 10]
    all_kmeans=[]
    for n_clusters in kmeans_clusters:
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
        y_kmeans = kmeans.predict(X)
        km_out = pd.Series(y_kmeans,index=X.index)
        km_out+=1
        km_out.name = 'K-Means Cluster (of %i)'%n_clusters
        
        all_kmeans.append(km_out)
    all_kmeans=pd.DataFrame(all_kmeans).transpose()
    # embed()
    new_analysis = pd.merge(new_analysis, all_kmeans, left_index=True, right_index=True, how='left')
    
    # embed()
    
    ''' KNN Additions '''
    if mda_file!='':
        X = zscores[zscores.columns[1:]]
        y = zscores['Group']
        y = y.replace(exp,1)
        y = y.replace(control,0)
        # embed()
        for n_clusters in kmeans_clusters:
            from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
            from sklearn.cluster import SpectralClustering
            # embed()
            # y = all_kmeans['K-Means Cluster (of %i)'%n_clusters]
            knn = KNeighborsClassifier(n_neighbors=n_clusters).fit(X,y)
            # nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
            # distances, indices = nbrs.kneighbors(X)
            embed()
            knn = pd.DataFrame(knn.predict(X), index=X.index)
            
            knn=knn.rename(columns={0:'kNN Cluster (of %i)'%n_clusters})
            
            new_analysis_tmp = pd.merge(new_analysis, knn, left_index=True, right_index=True, how='left')
            all_ranks = pd.DataFrame()
            tmp = new_analysis_tmp[['RFA MDA (5000 trees)','kNN Cluster (of %i)'%n_clusters]]
            
            # embed()
            
            for cluster in tmp['kNN Cluster (of %i)'%n_clusters].unique():
                tmp2 = tmp[tmp['kNN Cluster (of %i)'%n_clusters]==cluster]
                tmp2=tmp2.sort_values('RFA MDA (5000 trees)', ascending=False)
                tmp2['Rank']=range(1,len(tmp2)+1)
                all_ranks = all_ranks.append(tmp2)
            tmp2 = tmp.groupby('kNN Cluster (of %i)'%n_clusters).mean()
            all_ranks = all_ranks.rename(columns={'Rank':'MDA Rank in kNN Cluster (of %i)'%n_clusters})
            tmp2 = tmp2.rename(columns={'RFA MDA (5000 trees)':'Mean MDA/kNN Cluster (of %i)'%n_clusters})
            all_ranks = pd.merge(all_ranks,tmp2, left_on='kNN Cluster (of %i)'%n_clusters, right_index=True)
            
            # del all_ranks['kNN Cluster (of %i)'%n_clusters]
            del all_ranks['RFA MDA (5000 trees)']
            # print('STOPPED')
            # embed()
            new_analysis = pd.merge(new_analysis, all_ranks, left_index=True, right_index=True, how='left')
    # embed()
    # spc = SpectralClustering(n_clusters=10)
    # spc.fit(X)
    # res = spc.fit_predict(X)
    
   
    # embed()
    
    '''
    # MDA:
    # from: http://blog.datadive.net/selecting-good-features-part-iii-random-forests/
    
    X = samplesC.append(controlsC)
    Y = pd.Series([0]*len(X))
    Y.loc[samplesC.index]=1
    X = X.transpose().dropna().transpose()
    
    from sklearn.cross_validation import ShuffleSplit
    from sklearn.metrics import r2_score
    from collections import defaultdict
    from sklearn.ensemble import RandomForestRegressor
     
    rf = RandomForestRegressor(n_estimators=5000,random_state=123456, max_features=7)
    scores = defaultdict(list)
     
    #crossvalidate the scores on a number of different random splits of the data
    for train_idx, test_idx in ShuffleSplit(len(X), 10, .3):
        # embed()
        X_train, X_test = X.ix[train_idx], X.ix[test_idx]
        Y_train, Y_test = Y.ix[train_idx], Y.ix[test_idx]
        r = rf.fit(X_train, Y_train)
        acc = r2_score(Y_test, rf.predict(X_test))
        for c in X.columns:
            X_t = X_test.copy()
            tmp = X_t[c].values
            np.random.shuffle(tmp)
            X_t[c] = tmp
            shuff_acc = r2_score(Y_test, rf.predict(X_t))
            scores[c].append((acc-shuff_acc)/acc)
    print "Features sorted by their score:"
    print sorted([(round(np.mean(score), 4), feat) for
                  feat, score in scores.items()], reverse=True)
    embed()
    '''
    
    ''' ADD CORRELATIONS '''
    
    groups = data['Group']
    groups[groups==exp]=1
    groups[groups==control]=0
    from scipy.stats import spearmanr, pearsonr
    corr_out = []
    for c in data2.columns:
        # embed()
        if str(data2[c].mean())!='nan':
            vals = data2[c].dropna()
            group = groups.loc[vals.index,:]
            
            # adding kendall tau correlation (per Bob request 11/5/2020
            ktau, kpval = stats.kendalltau(group.values,vals.values) # note I am arbitrarily cutting off the last item... may not be fair
            
            corr_out.append({'MRM Name2':c,
                'Spearman_r': spearmanr(group.values, vals.values)[0], 
                'Spearman_p': spearmanr(group.values, vals.values)[1],
                'Pearson_r': pearsonr(group.values, vals.values)[0], 
                'Pearson_p': pearsonr(group.values, vals.values)[1],
                'Kendall_tau': ktau, 
                'Kendall_p': kpval
                
                
                
                # data.columns[1]+'_r': spearmanr(data[data.columns[1]].values, data[c].values)[0], 
                # data.columns[1]+'_p': spearmanr(data[data.columns[1]].values, data[c].values)[1], 
                # data.columns[2]+'_r': spearmanr(data[data.columns[2]].values, data[c].values)[0], 
                # data.columns[2]+'_p': spearmanr(data[data.columns[2]].values, data[c].values)[1], 
                
            })
            # embed()
        # else:
            # embed()
    corr_out = pd.DataFrame(corr_out)
    corr_out.index=corr_out['MRM Name2']
    
    
    
    
    
    ''' MERGING IT ALL TOGETHER '''  
    # embed()
    # out = pd.merge(table1_mets, out, left_on='MRM Name', right_on='MRM Name_x', how='outer')[['MRM Name','Pathway Name','VIP Score','P Value','FDR','q-Value','Power','MDA','Fold Change','Z Score']]
    # 12/12/2017: REMOVED MDA AND POWER:
    # out = pd.merge(table1_mets, out, left_on='MRM Name', right_on='MRM Name_x', how='outer')[['MRM Name','Pathway Name','VIP Score','P Value','FDR','q-Value','Fold Change','Z Score']]
    # 6/14/2018: removed linear Fold Change, calculated above - will be replaced with log lin FC calculated in new analysis, above
    
    ''' merge new analysis: (SEM, SD 1.96 ETC) '''
    # embed()
    out = pd.merge(out,new_analysis, left_on='MRM Name_y', right_index=True)
    
    out = pd.merge(out,logreg_out, left_on='index', right_on='LR MRM Name')
    
    # embed()
    ''' merge correlations '''
    out = pd.merge(out,corr_out, left_on='MRM Name_y', right_index=True)
    # embed()
    
    # adding columns per 9/24 request:
    out['Percent Difference (Case/Control -1)'] = out['Fold Change']-1
    out['Standard Error of the Difference (%)'] = out['Linear SEM']-1
    
    ''' merge it all together '''
    cols_to_include = [
        'MRM Name',
        'HMDB number',
        'CAS Number',
        'Pathway Name',
        'VIP Score',
        'Z Score',
        'Linear Z Score',
        'P Value',
        'Kruskal-Wallis p',
        # 'Mann-Whitney p (1-sided)', # don't included 1-sided per Bob Email 10/17/2018
        'Mann-Whitney U test p',
        'Welch\'s ttest p',
        'D\'Agostino & Pearson omnibus p',
        'D\'Agostino & Pearson omnibus p (Linear)',
        'Shapiro-Wilk normality p',
        'Shapiro-Wilk normality p (Linear)',
        'Control Outliers',
        'Control Low Outlier Threshold',
        'Control High Outlier Threshold',
        'Case Outliers',
        'Case Low Outlier Threshold',
        'Case High Outlier Threshold',
        'FDR',
        'q-Value',
        'Control Mean',
        'Control Median',
        'Control SD*',
        'Control SEM*',
        'Case Mean',
        'Case Median',
        'Case SD*',
        'Case SEM*',
        'Linear Control Mean (Geomean)',
        '2.50%','97.50%','Hi/Lo Ratio',
        'Linear Case Mean (Geomean)',
        'Fold Change',
        'Log2 SEM',
        'Linear SEM',
        'Percent Difference (Case/Control -1)',
        'Standard Error of the Difference (%)',
        'Spearman_r',
        'Spearman_p',
        'Pearson_r', 
        'Pearson_p',
        'Kendall_tau',
        'Kendall_p',
        'LogReg_AUC',
        'LogReg_AUC(split)',
        'LogReg_p',
        'LogReg_p(split)']
        
        
    if run_het==True:
        cols_to_include += [
            'Control Heteroscedasticity (linear space Breusch-Pagan p)',
            'Control Heteroscedasticity (log space Breusch-Pagan p)',
            'Case Heteroscedasticity (linear space Breusch-Pagan p)',
            'Case Heteroscedasticity (log space Breusch-Pagan p)'
        ]   
    
    
    
    if mda_file!='':
        cols_to_include.append('RFA MDA (5000 trees)')
        for n_clusters in kmeans_clusters:
            cols_to_include.append('kNN Cluster (of %i)'%n_clusters)
            cols_to_include.append('MDA Rank in kNN Cluster (of %i)'%n_clusters)
            cols_to_include.append('Mean MDA/kNN Cluster (of %i)'%n_clusters)
            cols_to_include.append('K-Means Cluster (of %i)'%n_clusters)
    
        
            
    if 'Mann-Whitney U test p' in table1_mets.columns:
        del table1_mets['Mann-Whitney U test p']
        del table1_mets['q-Value']
    # embed()
    out = out.rename(columns={'MRM Name':'dB MRM Name'})
    out = pd.merge(table1_mets, out, left_on='MRM Name', right_on='MRM Name_x', how='outer')
    # embed()
    del out['MRM Name']
    out = out.rename(columns={'MRM Name_x':'MRM Name'}) # this maps the original value of the MRM column to our output
    out = out[cols_to_include] # EDITED right merge index on 1/13/2021
    # out = pd.merge(table1_mets, out, left_on='MRM Name', right_on='index', how='outer')[cols_to_include] #
    # out = pd.merge(table1_mets, out, left_on='MRM Name', right_on='MRM Name', how='outer')[cols_to_include] # edited 2/9/21
    
    if len(out)>num_mets_measured:
        print('merge error, out >> measured metabolites')
        print('this still needs to be fixed - but its a problem with met mapping in mtdb')
        # first step is to remove duplicates
        # also look at len(table1_mets)
        embed()
    
    out.index = range(1,len(out)+1)

    import datetime as dt
    if VIP_cutoff=='':
        print('new code calls for listing VIP cutoff - see CFS2 run all')
        embed()
        
    # RENAME COLUNS PER BOB's request 9/24/2018    
    out = out.rename(columns={
        'VIP Score':'PLSDA VIP Score',
        'P Value':'ttest p',
        'Control Mean':'Log2 Control Mean',
        'Control Median':'Log2 Control Median',
        'Case Mean':'Log2 Case Mean',
        'Case Median':'Log2 Case Median',
        '2.50%':'2.50% Threshold#',
        '97.50%':'97.50% Threshold+'
    })
    if mda_file!='':
        out = out.rename(columns={
            'RFA MDA (5000 trees)':'Random Forest MDA (5000 trees)'
        })
        tmp = out['Random Forest MDA (5000 trees)'].sort_values(ascending=False)
        tmp = pd.DataFrame(tmp)
        tmp['Random Forest Rank'] = range(1,len(tmp)+1)
        del out['Random Forest MDA (5000 trees)']
        out = pd.merge(out, tmp, left_index=True, right_index=True)
        
    
    # renaming cols per Bob's request - 2-11-2019
    out = out.rename(columns={
        'Z Score':'Log Space Z-Score',
        'ttest p':'Student\'s ttest p',
        # 'kNN Cluster':'kNN Cluster (of 10)'
    })
    
    # embed()
    
    
    def package_mdas(knn_mda_in, n_clusters):
        # embed()
        knn_mda2 = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).sum()
        knn_mda2['Count per cluster'] = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']
        positive_mdas = knn_mda_in[knn_mda_in['Random Forest MDA (5000 trees)']>0]
        # embed()
        knn_mda2['Count of positive MDAs per cluster'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']
        knn_mda2['Mean Z-score'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).mean()['Log Space Z-Score']
        knn_mda2['SD'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).std()['Log Space Z-Score']
        # embed()
        knn_mda2['Sum of positive MDAs per cluster'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).sum()['Random Forest MDA (5000 trees)']
        knn_mda2['Fraction of positive']=knn_mda2['Sum of positive MDAs per cluster']/knn_mda2['Sum of positive MDAs per cluster'].sum()
        
        knn_mda2 = knn_mda2.sort_values('Fraction of positive', ascending=False)
        knn_mda2['Rank']=range(1,len(knn_mda2)+1)
        
        for i in range(1,11):
            cluster_mdas = positive_mdas[positive_mdas['kNN Cluster (of %i)'%n_clusters]==i]
            cluster_mdas = cluster_mdas.sort_values('Random Forest MDA (5000 trees)', ascending=False)[:20]
            # embed()
            count=1
            
            for ii in cluster_mdas.index:
                name = cluster_mdas.loc[ii,'MRM Name']
                knn_mda2.loc[i,count] = name
                count+=1
        knn_mda2 = knn_mda2[[u'Random Forest MDA (5000 trees)', u'Log Space Z-Score',
           u'Count per cluster', u'Count of positive MDAs per cluster',
           u'Sum of positive MDAs per cluster',
           u'Fraction of positive', u'Mean Z-score', u'SD', u'Rank']+range(1,21)]
        return knn_mda2
        
    def package_knn_vip(knn_mda_in, n_clusters):
        # embed()
        knn_mda2 = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).sum()
        knn_mda2['Count per cluster'] = knn_mda_in.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']

        positive_mdas = knn_mda_in[knn_mda_in['PLSDA VIP Score']>knn_vip_threshold]
        # embed()
        knn_mda2['Mean Z-score'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).mean()['Log Space Z-Score']
        knn_mda2['SD'] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).std()['Log Space Z-Score']
        # positive_mdas = knn_mda_in
        knn_mda2['Count of VIP>%0.1f per cluster'%(knn_vip_threshold)] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).count()['MRM Name']
        # embed()
        knn_mda2['Sum of VIP>%0.1f per cluster'%(knn_vip_threshold)] = positive_mdas.groupby('kNN Cluster (of %i)'%n_clusters).sum()['PLSDA VIP Score']

        knn_mda2['Fraction of positive']=knn_mda2['Sum of VIP>%0.1f per cluster'%(knn_vip_threshold)]/knn_mda2['Sum of VIP>%0.1f per cluster'%(knn_vip_threshold)].sum()

        knn_mda2 = knn_mda2.sort_values('Fraction of positive', ascending=False)
        knn_mda2['Rank']=range(1,len(knn_mda2)+1)
                
        for i in range(1,11):
            cluster_mdas = positive_mdas[positive_mdas['kNN Cluster (of %i)'%n_clusters]==i]
            cluster_mdas = cluster_mdas.sort_values('PLSDA VIP Score', ascending=False)[:20]
            # embed()
            count=1
            for ii in cluster_mdas.index:
                name = cluster_mdas.loc[ii,'MRM Name']
                knn_mda2.loc[i,count] = name
                count+=1
        del knn_mda2['Log Space Z-Score']
        # embed()
        num_clusters = len(knn_mda2.columns)-8
        knn_mda2 = knn_mda2[[u'PLSDA VIP Score', u'Count per cluster',
           u'Count of VIP>%0.1f per cluster'%(knn_vip_threshold), u'Sum of VIP>%0.1f per cluster'%(knn_vip_threshold),
           u'Fraction of positive', 'Mean Z-score','SD',u'Rank']+range(1,num_clusters+1)]
     
        return knn_mda2
        
    
    if mda_file!='':
        knn_mda_files = {}
        knn_mda_vip_files = {}
        for n_clusters in kmeans_clusters:
            knn_outfile = outfile.replace('TableS1','knnClusters_%i'%n_clusters)
            writer2 = pd.ExcelWriter(knn_outfile)
            
            # embed()
            
            
            knn_mda = out[['MRM Name','Random Forest MDA (5000 trees)','kNN Cluster (of %i)'%n_clusters,'Log Space Z-Score']]
            
            # addition 7/23/2020 to made a new table for same analysis but with vip>1
            knn_mda_vip = out[['MRM Name','kNN Cluster (of %i)'%n_clusters,'PLSDA VIP Score','Log Space Z-Score']]# filter for those with vip>1
            # knn_mda_vip = knn_mda_vip[knn_mda_vip['PLSDA VIP Score']>1]
            # embed()
            knn_mda_vip_packaged = package_knn_vip(knn_mda_vip, n_clusters)
            # end new analysis 7/23/2020
                
            knn_mda_packaged = package_mdas(knn_mda, n_clusters)
            knn_mda_files[str(n_clusters)] = knn_mda_packaged
            knn_mda_vip_files[str(n_clusters)] = knn_mda_vip_packaged
            
            out['kNN Cluster (of %i)'%n_clusters] = out['kNN Cluster (of %i)'%n_clusters].fillna(99)
            # embed()
                
            knn_groups = out['kNN Cluster (of %i)'%n_clusters].unique().tolist()
            knn_groups.sort()
            # embed()
            for knn_group in knn_groups:
                groupmets = out[out['kNN Cluster (of %i)'%n_clusters]==knn_group]
                # embed()
                groupmets = groupmets[['MRM Name','PLSDA VIP Score','Log Space Z-Score','Mann-Whitney U test p','FDR','Random Forest MDA (5000 trees)','Pathway Name','kNN Cluster (of %i)'%n_clusters]]
                # try:
                # embed()
                groupmets.to_excel(writer2, str(int(knn_group)))
                # except:
                    # embed()
            # knn_mda = knn_mda2
            
            ''' OLD '''
            '''
            
            knn_groups = out['kNN Cluster'].unique().tolist()
            knn_groups.sort()
            for knn_group in knn_groups:
                groupmets = knn[knn['kNN Cluster']==knn_group]
                groupmets['new_index'] = groupmets.index.str.replace('(','').str.replace(')','').str.replace(':','')  
                groupmets = pd.merge(groupmets, mda, left_on='new_index', right_index=True)
                groupmets = pd.merge(groupmets, table1_mets, left_index=True, right_on='MRM Name')
                groupmets = groupmets[['MRM Name','RFA MDA (5000 trees)','Pathway Name','kNN Cluster']]
                groupmets.to_excel(writer2, str(knn_group))
            '''
            
            
            ''' kMeans groups '''
            '''
            km_groups = km_out.unique().tolist()
            km_groups.sort()
            # embed()
            for km_group in km_groups:
                groupmets = km_out[km_out==km_group]
                groupmets = pd.DataFrame(groupmets)
                groupmets['new_index'] = groupmets.index.str.replace('(','').str.replace(')','').str.replace(':','')
                # embed()
                groupmets = pd.merge(groupmets, mda, left_on='new_index', right_index=True)
                groupmets = pd.merge(groupmets, table1_mets, left_index=True, right_on='MRM Name')
                groupmets = groupmets[['MRM Name','RFA MDA (5000 trees)','Pathway Name','K-Means Cluster (of 10)']]
                groupmets.to_excel(writer2, str(km_group))
            '''
            
            writer2.save()
        
        # embed()
        
    
    
    writer = pd.ExcelWriter(outfile)
    readme = {'Study':study_name,
    'Group':exp,
    'Database File':db_file,
    'Analysis date':dt.datetime.today().strftime("%m/%d/%Y"),
    'VIP Score Cutoff':VIP_cutoff,
    'N (control=%s)'%control:len(controls),
    'N (case=%s)'%exp:len(samples),
    'Metabolites Measured':num_mets_measured,
    'Metabolites Targeted':num_mets_targeted,
    'Methods':'Outlier calculation was by the method of Pierce implemented in Python',
    'Methods ':'Log2 SEM calculation was by dividing the standard deviation of the original distribution by the square root of the sample size for Log2 transformed data',
    'Methods  ':'Linear SEM calculation was by dividing the standard deviation of the original distribution by the square root of the sample size for Linear data',
    'Notes':'All p values are for 2-sided tests unless otherwise specified',
    # 'ML iterations':ML_iter,
    }
    readme = pd.DataFrame(pd.Series(readme),columns=['Table S1'])
    readme.to_excel(writer,'README')
    zscores.to_excel(writer,'Zscores')
    out.to_excel(writer,'TableS1')
    if mda_file!='':
        for n_clusters in kmeans_clusters:
            knn_mda_files[str(n_clusters)].to_excel(writer,'kNN (%s) Ranking by MDA'%str(n_clusters))
            knn_mda_vip_files[str(n_clusters)].to_excel(writer,'kNN (%s) Ranking by VIP>%0.1f'%(n_clusters,knn_vip_threshold))
    writer.save()
    # embed()
    
#!/usr/bin/python
#
# peirce_dev.py
# created 16 Jul 2013
# updated 23 Oct 2014
#
#### MODULES ####
import numpy
import scipy.special




#### FUNCTION ####


def get_zscores(data, case, control):
    # Log transform data
    from scipy import stats
    # embed()
    log2 = np.log2(data[data.columns[1:]])
    log2['Group']=data['Group']
    data = log2

    count = data.count()
    full_rank = count[count==count.max()]
    other = count[(count!=0)&(count!=count.max())]
    other = data[other.index]
    other['Group'] = data['Group']
    other = other.dropna()
    data = data[full_rank.index]
    data['Group'] = data['Group']
    
    def get_z(data, case, control):
        
        
        # split case and controls
        control = data[data['Group']==control][data.columns[:-1]]
        case = data[data['Group']==case][data.columns[:-1]]

        # calculate z-scores
        control_z = pd.DataFrame(stats.zscore(control, ddof=1), index=control.index, columns=control.columns)
        control_mean = control.mean()
        control_std = control.std()
        case_z = (case-control_mean)/control_std
        return control_z.append(case_z)
        
    full = get_z(data, case, control)
    other = get_z(other, case, control)
    # embed()
    # data = full['case'].append(full['control'])
    
    out = pd.merge(full, other, left_index=True, right_index=True, how='outer')
    out['Group'] = data['Group']
    case_z = out[out['Group']==case]
    control_z = out[out['Group']==control]
    
    return {'case':case_z, 'control':control_z}

def peirce_dev(N, n, m):
   """
   Name:     peirce_dev
   Input:    - int, total number of observations (N)
             - int, number of outliers to be removed (n)
             - int, number of model unknowns (m)
   Output:   float, squared error threshold (x2)
   Features: Returns the squared threshold error deviation for outlier 
             identification using Peirce's criterion based on Gould's
             methodology
   """
   # Assign floats to input variables:
   N = float(N)
   n = float(n)
   m = float(m)
   #
   # Check number of observations:
   if N > 1:
      # Calculate Q (Nth root of Gould's equation B):
      Q = (n**(n/N)*(N - n)**((N - n)/N))/N
      #
      # Initialize R values (as floats)
      Rnew = 1.0  
      Rold = 0.0  # <- Necessary to prompt while loop
      #
      # Start iteration to converge on R:
      while ( abs(Rnew - Rold) > (N*2.0e-16) ):
         # Calculate Lamda 
         # (1/(N-n)th root of Gould's equation A'):
         ldiv = Rnew**n
         if ldiv == 0:
            ldiv = 1.0e-6
         Lamda = ((Q**N)/(ldiv))**(1.0/(N - n))
         #
         # Calculate x-squared (Gould's equation C):
         x2 = 1.0 + (N - m - n)/n*(1.0 - Lamda**2.0)
         #
         # If x2 goes negative, return 0:
         if x2 < 0:
            x2 = 0.0
            Rold = Rnew
         else:
            # Use x-squared to update R (Gould's equation D):
            Rold = Rnew
            Rnew = (
               numpy.exp((x2 - 1)/2.0)*
               scipy.special.erfc(numpy.sqrt(x2)/numpy.sqrt(2.0))
               )
         #
   else:
      x2 = 0.0
   return x2
   
def ts1_summary_fig(ts1_file, outdir='', sort_col='PLSDA VIP Score'):
    # NOTE this could be a good example to follow: https://matplotlib.org/3.1.1/gallery/misc/table_demo.html#sphx-glr-gallery-misc-table-demo-py
    import pylab as plt
    import seaborn as sns
    data = pd.read_excel(ts1_file, 'TableS1')
    
    data['Pathway Name'] = data['Pathway Name'].str.replace(' Metabolism','s')
    data['Pathway Name'] = data['Pathway Name'].str.replace("\)s",")")
    
    if sort_col=='PLSDA VIP Score':
        data = data.sort_values('PLSDA VIP Score', ascending=True)
        type='VIP'
    elif sort_col=='Random Forest MDA (5000 trees)':
        data = data.sort_values('Random Forest Rank', ascending=False)
        type='MDA'
    data = data[-30:]
    sorted_pwy = data.groupby('Pathway Name').count()['MRM Name'].sort_values(ascending=False)
    sorted_pwy.name = 'Pathway Rank'
    sorted_pwy = pd.DataFrame(sorted_pwy)
    print(sorted_pwy)
    
    # embed()
    
    data['Metabolite Rank'] = range(1,31)
    data = data[['MRM Name', 'Metabolite Rank','Pathway Name','PLSDA VIP Score','Random Forest MDA (5000 trees)','Random Forest Rank']]
    data = pd.merge(data, sorted_pwy, left_on='Pathway Name', right_index=True)
    data = data.sort_values('Pathway Rank', ascending=False)
    # embed()
    
    '''
    unique_pathways = data['Pathway Name'].unique()
    # colors = sns.color_palette("hls", len(unique_pathways))
    colors = sns.color_palette("Wistia", len(unique_pathways))
    # colors = sns.color_palette("Paired")
    # plt.show(colors)
    # colors = sns.color_palette("Blues", len(unique_pathways))
    # embed()
    colors = dict(zip(unique_pathways,colors))
    data['color'] = data['Pathway Name'].map(colors)
    '''
    
    # embed()
    # data = data.sort_values('Pathway Name')
    
    # data['size']=100
    sns.set(rc={'figure.figsize':(10,8)})
    # embed()
    
    # ax = sns.scatterplot(data=data, x=sort_col, y="Metabolite Rank", hue="Pathway Name", s=100, palette='bright')
    palette = sns.color_palette("bright")
    palette+=[(0,0,0),(1,1,1),(0.2,0.2,0.2),(0.4,0.4,0.4),(0.6, 0.6, 0.6),(0.8,0.8,0.8)]
    if len(sorted_pwy)>len(palette):
        print('more pathways than colors')
        embed()
    # embed()
    ax = sns.scatterplot(data=data, x=sort_col, y="Metabolite Rank", hue="Pathway Name", s=100, palette=palette[:len(sorted_pwy)])
    # ax = sns.scatterplot(data=data, x=sort_col, y="Metabolite Rank", hue="Pathway Name", s=100, palette=sns.color_palette("hls", len(sorted_pwy)))
    # ax = sns.scatterplot(data=data, x=sort_col, y="Metabolite Rank", hue="Pathway Name", s=100, palette=sns.color_palette("Paired")[:len(sorted_pwy)])
    # plt.setp(ax.get_legend().get_texts(), fontsize='22') # for legend text
    # plt.setp(ax.get_legend().get_title(), fontsize='32') # for legend title
    # ax.set_xticklabels(ax,res.get_xmajorticklabels(), fontsize = 18)
    
    # embed()
    # line=0
    # for a in ax.get_legend().get_texts():
        # if line!=0:
            # a.set_text(sorted_pwy.index[line-1])
        # line+=1
        
    l = ax.get_xlabel()
    ax.set_xlabel(l, fontsize=20)
    l = ax.get_ylabel()
    ax.set_ylabel(l, fontsize=20)
    # embed()
    data = data.sort_values('Metabolite Rank')
    plt.yticks(range(1,31),data['MRM Name'].tolist())
    plt.tight_layout()
    # plt.show()
    
    
    filename = ts1_file.split('\\')[-1].split('/')[-1]
    # embed()
    plt.title(outdir + type + ' ' + filename)
    plt.savefig(outdir + type + '_' + filename.replace('.xlsx','.png'))
    plt.savefig(outdir + type + '_' + filename.replace('.xlsx','.pdf'))
    plt.savefig(outdir + type + '_' + filename.replace('.xlsx','.svg'))
    plt.close()
    
    '''
    
    plt.figure(figsize=(5,7))
    points = plt.scatter(data[sort_col], data['y_coord'], c=data['color'], s=100, label=data['Pathway Name'].tolist())
    plt.xlabel(sort_col)
    plt.yticks(range(1,31),data['MRM Name'].tolist())
    plt.tight_layout()
    plt.legend([points], [data['Pathway Name'].tolist()])
    # plt.legend(handles=[points])
    plt.show()
    embed()
    '''
   
   
def compare_clusters():
    # check that my kmeans clusters are similar to the ones for metaboanalys
    mine = pd.read_excel(r'D:\Dropbox (Personal)\naviauxlab_informatics\GWI\Results/GWI_new_GWI_TableS1_mtdb4.4.2.xlsx','TableS1')
    mine = mine[['MRM Name',u'K-Means Cluster (5)']]
    metabo = pd.read_csv(r'D:\Dropbox (Personal)\naviauxlab_informatics\GWI\metaboanalyst/metaboanalyst_roc_univ.csv',index_col=0)
    metabo = metabo['clusters'].reset_index()
    # metabo = metabo.rename(columns={'index','MRM Name'})
    mine['index'] = mine['MRM Name'].str.replace(':','')
    mine['index'] = mine['index'].str.replace('(','')
    mine['index'] = mine['index'].str.replace(')','')
    
    out = pd.merge(mine,metabo)
    out.to_csv('kmean_test.csv')
    
    embed()
    
def check_pval_qval():
    # check to see why FDR and qval are exact in study 5 but not 1 or 2
    data1 = pd.read_excel(r'D:\Dropbox (Personal)\naviauxlab_informatics\Dutch Depression\RMDD\1\Results\1_Male_RMDD-IR_vs_Controls_RMDD-IR_TableS1_mtdb5.0.2.xlsx','TableS1')
    data2 = pd.read_excel(r'D:\Dropbox (Personal)\naviauxlab_informatics\Dutch Depression\RMDD\5\Results\5_Controls_Males_v_Females_Male Controls_TableS1_mtdb5.0.2.xlsx','TableS1')
    
    pv1 = data1['P Value']
    pv2 = data2['P Value']
    
    import numpy as np 
    from . import qvalue

    # pv = np.random.uniform(0.0, 1.0, size = (1000,))
    pv = pd.Series(np.random.normal(0.5, .1, size = (1000,)))
    # embed()
    
    def get_pqf(pv):
        pv = pv.sort_values()
        pv.index=range(0,len(pv))
        qv = qvalue.estimate(pv)
        
        
        from statsmodels.stats.multitest import fdrcorrection_twostage

        res = fdrcorrection_twostage(pv, alpha=0.05, method='bh')
        rej = res[0]
        pval_corr = res[1]
        
        
        
        import pylab as plt
        # embed()
        
        plt.scatter(pv,pv.index)
        plt.show()
        
        pd.Series(pv).hist(color='red', bins=25)
        plt.xlim([0,1])
        plt.show()
        pd.Series(qv).hist(color='blue', bins=25)
        plt.xlim([0,1])
        plt.show()
        pd.Series(pval_corr).hist(color='green', bins=25)
        plt.xlim([0,1])
        plt.show()
    
        tmp = pd.DataFrame([pv,qv,pval_corr], index=['pvals','qvals','fdr']).transpose()
        return tmp
    # test = get_pqf(pv)
    s1 = get_pqf(pv1)
    s2 = get_pqf(pv2)
    
    embed()
   
if __name__=='__main__':
    # compare_clusters()
    # check_pval_qval()
    main()
    # pass