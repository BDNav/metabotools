import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from IPython import embed

from metabotools.metabotools import load_db, series2dataframe

# table1 = pd.read_excel('Females/DEPRESS_tableOne_VIP_Score1.2.xlsx','Metabolites')
# table1.to_csv('Females/figures/diag.csv')


data = ''
group = ''
control=''
sample=''
Z_score = ''
db=''
alias_db=''
global_diag = ''
VIP_score = ''


def init(full_data, control_l, sample_l, Z_score_l, db_file, table1_mets, VIP_score_l):
    global data, group, control, sample, Z_score, z, alias_db, db, global_diag, VIP_score
    control = control_l
    sample = sample_l
    Z_score = Z_score_l
    full_data = pd.read_csv(full_data).dropna(axis=1, how='all')
    data = full_data[~full_data['Sample Name'].str.contains('REP')].set_index('Sample Name').drop('Group', axis=1)
    group = full_data[~full_data['Sample Name'].str.contains('REP')].set_index('Sample Name')['Group']
    db, alias_db = load_db(db_file)
    global_diag = table1_mets
    VIP_score = VIP_score_l
    
def ctrl(data, group):
    return data[group == control]

def exp(data, group):
    return data[group == sample]

def effect_size(data, mean, std):
    return (data - mean) / std

def hilo_plot(outfile, show=0):
    from scipy.stats import ttest_ind
    ''' was trying to see why hilo doesn't line up with prev figure - should just remove it...
    z_control = effect_size(np.log2(data), np.log2(ctrl(data, group)).mean(), np.log2(ctrl(data, group)).std())
    z_exp = effect_size(np.log2(data), np.log2(exp(data, group)).mean(), np.log2(exp(data, group)).std())
    
    z_control = z_control.loc[group[group==control].index]
    z_exp = z_exp.loc[group[group==sample].index]
    z = z_control.append(z_exp)
    # embed()
    '''
    z = effect_size(np.log2(data), np.log2(ctrl(data, group)).mean(), np.log2(ctrl(data, group)).std())
    
    plt_data = pd.DataFrame({'hi': (z >= Z_score).sum(axis=1), 'lo': (z <= -1*Z_score).sum(axis=1)})
    plt_data['total'] = plt_data.hi + plt_data.lo
    
    plt_data = pd.concat([plt_data, group], axis=1)
    plt_sample = plt_data[plt_data.Group==sample]
    plt_control = plt_data[plt_data.Group==control]
    plt_sample = plt_sample.sort_values('total', ascending=False)
    plt_control = plt_control.sort_values('total', ascending=False)
    plt_data = plt_sample.append(plt_control)
    # .sort_values(['Group', 'total'], ascending=[True, False])
    # '''
    # embed()
    fig, ax = plt.subplots(figsize=(20,15))
    fontsize = 30
    width = 0.35
    bar_positions = np.arange(len(plt_data)) + .35

    ax.set_xlim(0, max(bar_positions) + 1)

    ax.bar(bar_positions, plt_data.lo, color='g', width=width, label='Low')
    ax.bar(bar_positions, plt_data.hi, color='r', bottom=plt_data.lo, width=width, label='High')
    ax.axhline(y=plt_data[plt_data.Group == control]['total'].median(), color='b', linestyle='--', linewidth=3, label='%s Median'%control)
    ax.axhline(y=plt_data[plt_data.Group == sample]['total'].median(), color='y', linestyle='--', linewidth=3, label='%s Median'%sample)
    
    control_count = len(plt_data[plt_data.Group == control]['total'])
    sample_count = len(plt_data[plt_data.Group == sample]['total'])
    
    control_mean = plt_data[plt_data.Group == control]['total'].mean()
    sample_mean = plt_data[plt_data.Group == sample]['total'].mean()
    
    control_sem = plt_data[plt_data.Group == control]['total'].sem()
    sample_sem = plt_data[plt_data.Group == sample]['total'].sem()
    
    pval_total = ttest_ind(plt_data[plt_data.Group == sample]['total'], plt_data[plt_data.Group == control]['total'])
    pval_lo = ttest_ind(plt_data[plt_data.Group == sample]['lo'], plt_data[plt_data.Group == control]['lo'])
    
    control_lo_mean = plt_data[plt_data.Group == control]['lo'].mean()
    sample_lo_mean = plt_data[plt_data.Group == sample]['lo'].mean()
    
    control_lo_sem = plt_data[plt_data.Group == control]['lo'].sem()
    sample_lo_sem = plt_data[plt_data.Group == sample]['lo'].sem()
    
    if pval_total[1]<0.05:
        text = "%s Total: Mean +/- SEM = %0.0f +/- %0.1f; pval < %0.0E"%(sample,sample_mean, sample_sem, pval_total[1])
    else:
        text = "%s Total: Mean +/- SEM = %0.0f +/- %0.1f; pval: ns"%(sample,sample_mean, sample_sem)
    if pval_lo[1]<0.05:
        text+= "\n%s Low: Mean +/- SEM = %0.0f +/- %0.1f; pval < %0.0E"%(sample,sample_lo_mean, sample_lo_sem, pval_lo[1])
    else:
        text+= "\n%s Low: Mean +/- SEM = %0.0f +/- %0.1f; pval: ns"%(sample,sample_lo_mean, sample_lo_sem)
    plt.text(2, plt_data[plt_data.Group == sample]['total'].max()+3,text, fontsize=fontsize*.8)
    plt_data.to_csv(outfile.split('.')[0]+'.csv')
    print(text)
    
    text =  "%s Total: Mean +/- SEM = %0.0f +/- %0.1f"%(control,control_mean, control_sem)
    text+= "\n%s Low: Mean +/- SEM = %0.0f +/- %0.1f"%(control,control_lo_mean, control_lo_sem)
    plt.text(len(plt_data[plt_data.Group == sample]['total']), plt_data[plt_data.Group == control]['total'].max()+3,text, fontsize=fontsize*.8)
    
    ax.set_xticks(bar_positions + width/2)
    ax.set_xticklabels(plt_data.index.tolist(), rotation=90, fontsize=fontsize*0.7)
    ax_ticks = range(0,plt_data.max()['total']+10,10)
    ax.set_yticks(ax_ticks)
    ax.set_yticklabels(ax_ticks, fontsize=fontsize)
    ax.set_ylabel('Metabolites out of normal range', fontsize=fontsize)
    ax.set_title('Hi-Lo Metabolite Plot Z=%0.2f (%s n = %i; %s n=%i)'%(Z_score,sample,sample_count,control,control_count), fontsize=fontsize)
    ax.tick_params(bottom="off", top="off", left="off", right="off")
    [v.set_visible(False) for k,v in ax.spines.items()]

    plt.legend(fontsize=fontsize*.8)
    if show==0:
        plt.savefig(outfile)    
        plt.close()
    else:
        plt.show()


def plot_horz_pathway(dir, s, db2, name):
    print(name)
    z = pd.merge(db2, series2dataframe(s, index_name='MRM Index', values_name='data'), on='MRM Index')
    z['hi'] = z['data'] >= Z_score
    z['lo'] = z['data'] <= -1* Z_score

    fig_order = z.groupby('Pathway Name').max()['Figure Order']
    plt_data = z.groupby('Pathway Name').sum().drop('Figure Order', axis=1)
    plt_data = pd.concat([plt_data, fig_order], axis=1).sort_values('Figure Order', ascending=False)

    fig, ax = plt.subplots(figsize=(10,7.5))
    plt.subplots_adjust(left=.3, right=.95)

    width = 0.35
    bar_positions = np.arange(len(plt_data)) + .35

    ax.set_ylim(0, max(bar_positions) + 1)

    ax.barh(bar_positions, plt_data.lo, color='g', height=width, label='Low')
    ax.barh(bar_positions, plt_data.hi, color='r', left=plt_data.lo, height=width, label='High')

    ax.set_yticks(bar_positions + width/2)
    ax.set_yticklabels(plt_data.index.tolist(), fontsize=6, va='center')

    ax.set_xticks([0,20,40])
    ax.set_xlabel('Metabolites out of normal range')
    ax.set_title('Hi-Lo Pathway Plot - {}'.format(name))
    ax.tick_params(bottom="off", top="off", left="off", right="off")
    #[v.set_visible(False) for k,v in ax.spines.items()]
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.legend()
    plt.savefig(dir+'hilo_path_%s_z=%0.2f.png'%(name,Z_score))
    plt.close()
    

def plot_all_indiv_horz_pathways(dir):
    z = effect_size(np.log2(data), np.log2(ctrl(data, group)).mean(), np.log2(ctrl(data, group)).std())
    z = z.rename(columns=dict(zip(alias_db['MRM Name'], alias_db['MRM Index'])))
    db2 = db[db['MRM Index'].isin(z.columns.tolist())]
    for i in z.index.tolist():
        plot_horz_pathway(dir, z.loc[i], db2, i)
    

    

def plot_all_diag(outfile, DEBUG=False, yaxis=range(0,100,10)):
    z = effect_size(np.log2(data), np.log2(ctrl(data, group)).mean(), np.log2(ctrl(data, group)).std())
    z = z.rename(columns=dict(zip(alias_db['MRM Name'], alias_db['MRM Index'])))
    
    plt_data = pd.DataFrame({'hi': (z >= Z_score).sum(axis=1), 'lo': (z <= -1*Z_score).sum(axis=1)})
    plt_data['total'] = plt_data.hi + plt_data.lo
    plt_data = pd.concat([plt_data, group], axis=1)
    
    # plt_data = pd.concat([plt_data, group], axis=1)
    plt_sample = plt_data[plt_data.Group==sample]
    plt_control = plt_data[plt_data.Group==control]
    plt_sample = plt_sample.sort_values('total')
    plt_control = plt_control.sort_values('total')
    plt_data = plt_sample.append(plt_control)
    # embed()
    
    
    # .sort_values(['Group', 'total'], ascending=[True, False])
        
    control_index = plt_data[plt_data.Group==control].index 
    sample_index = plt_data[plt_data.Group==sample].index
    z_w_group = z.copy()
    z_w_group.loc[control_index,'Group'] = control
    z_w_group.loc[sample_index,'Group'] = sample
    control_mean = z_w_group[z_w_group.Group==control].mean()
    sample_mean = z_w_group[z_w_group.Group==sample].mean()
    # embed()
    def select_diag_mets(s, DEBUG=False):

        d = global_diag
        # embed()
        d['MRM Index'] = d['MRM Name'].map(dict(zip(alias_db['MRM Name'], alias_db['MRM Index'])))
        d = d.set_index('MRM Index')
        # embed()
        # d = pd.merge(d, pd.DataFrame(s), left_on='MRM Name', right_index=True)
        d = pd.merge(d, pd.DataFrame(s), left_index=True, right_index=True)
        d['abs'] = abs(d[s.name])
        
        total=d[d['abs']>Z_score]
        # embed()
       
        potential_diag = total[total['VIP Score']>VIP_score]
        print(s.name)
        print('total (val > z_score):', len(total))
        print('potential diag (val > VIP score):', len(potential_diag))
        print()
       
        
        if s.name in control_index:
            type='control'
        elif s.name in sample_index:
            type='sample'
        diag = []
        # embed()
        for i in potential_diag.index:
            # embed()
            # name = potential_diag.loc[i,'MRM Name']
            name = i #change me based on format of s
            val = potential_diag.loc[i,s.name]
                        
            if DEBUG==True:
            
                diag_true=False
                plt.boxplot([z_w_group[z_w_group.Group==control][i],z_w_group[z_w_group.Group==sample][i]])
                plt.scatter([1]*len(z_w_group[z_w_group.Group==control][i]),z_w_group[z_w_group.Group==control][i])
                plt.scatter([2]*len(z_w_group[z_w_group.Group==sample][i]),z_w_group[z_w_group.Group==sample][i])
                if type=='control':
                    plt.scatter(1,val,color='red')
                elif type=='sample':
                    plt.scatter(2,val,color='red')
            
            print(type, name, val, control_mean.loc[name], sample_mean.loc[name], )
            
            if type=='control': 
                # if control_mean.loc[name]<0 and val<0:
                if control_mean.loc[name]<sample_mean.loc[name] and val<0:
                    diag.append(name)
                    print ('diag')
                    diag_true=True
                # elif control_mean.loc[name]>0 and val>0:
                elif control_mean.loc[name]>sample_mean.loc[name] and val>0:
                    diag.append(name)
                    print ('diag')
                    diag_true=True
            elif type=='sample':
                # if sample_mean.loc[name]<0 and val<0:
                if sample_mean.loc[name]<control_mean.loc[name] and val<0:
                    diag.append(name)
                    print ('diag')
                    diag_true=True
                # elif sample_mean.loc[name]>0 and val>0:
                elif sample_mean.loc[name]>control_mean.loc[name] and val>0:
                    diag.append(name)
                    print ('diag')
                    diag_true=True
                    
            if DEBUG == True:
                if diag_true==True:
                    plt.title('DIAGNOSTIC ' + s.name + ': ' + i)
                else:
                    plt.title('NOT diag ' + s.name + ': ' + i)
                plt.show()
                plt.close()
            print ()
        
        
            
        print(s.name, 'diag:', len(diag), 'individual:', len(total)-len(diag))
        # embed()
        return {'name':s.name, 'diag':len(diag), 'ind':len(total)-len(diag)}
    
    diag_mets = [select_diag_mets(z.loc[x], DEBUG=DEBUG) for x in z.index.tolist()]

    plt_data = pd.DataFrame(diag_mets).set_index('name')
    plt_data['total'] = plt_data.diag + plt_data.ind
    plt_data = pd.concat([plt_data, group], axis=1)
    
    plt_sample = plt_data[plt_data.Group==sample]
    plt_control = plt_data[plt_data.Group==control]
    plt_sample = plt_sample.sort_values('total', ascending=False)
    plt_control = plt_control.sort_values('total', ascending=False)
    plt_data = plt_sample.append(plt_control)
    
    # .sort_values(['Group', 'total'], ascending=[True, False])
    
    ctmean = plt_data[plt_data.Group==control]['total'].mean()
    ctsem = plt_data[plt_data.Group==control]['total'].sem()

    print('control: mean total mets', ctmean, '+/-', ctsem)

    cdmean = plt_data[plt_data.Group==control]['diag'].mean()
    cdsem = plt_data[plt_data.Group==control]['diag'].sem()

    print('control: mean diag mets', cdmean, '+/-', cdsem)

    cimean = plt_data[plt_data.Group==control]['ind'].mean()
    cisem = plt_data[plt_data.Group==control]['ind'].sem()

    print('control: mean ind mets', cimean, '+/-', cisem)

    stmean = plt_data[plt_data.Group==sample]['total'].mean()
    stsem = plt_data[plt_data.Group==sample]['total'].sem()
    
    sdmean = plt_data[plt_data.Group==sample]['diag'].mean()
    sdsem = plt_data[plt_data.Group==sample]['diag'].sem()

    from scipy.stats import ttest_ind
    from scipy import stats

    # total_p = ttest_ind(plt_data[plt_data.Group==sample]['total'], plt_data[plt_data.Group==control]['total'])
    total_p = stats.mannwhitneyu(plt_data[plt_data.Group==sample]['total'], plt_data[plt_data.Group==control]['total'], use_continuity=True, alternative='two-sided') 
    # total_p = stats.mannwhitneyu(plt_data[plt_data.Group==sample]['total'], plt_data[plt_data.Group==control]['total']) 
    
    # diag_p = ttest_ind(plt_data[plt_data.Group==sample]['diag'], plt_data[plt_data.Group==control]['diag'])
    
    diag_p = stats.mannwhitneyu(plt_data[plt_data.Group==sample]['diag'], plt_data[plt_data.Group==control]['diag'], use_continuity=True, alternative='two-sided') 
    # diag_p = stats.mannwhitneyu(plt_data[plt_data.Group==sample]['diag'], plt_data[plt_data.Group==control]['diag']) 
        
    # embed()

    print('sample: mean diag mets', sdmean, '+/-', sdsem, 'p-val', diag_p[1])



    simean = plt_data[plt_data.Group==sample]['ind'].mean()
    sisem = plt_data[plt_data.Group==sample]['ind'].sem()
    
    control_count = len(plt_data[plt_data.Group==control])
    sample_count = len(plt_data[plt_data.Group==sample])

    # ind_p = ttest_ind(plt_data[plt_data.Group==sample]['ind'], plt_data[plt_data.Group==control]['ind'])
    ind_p = stats.mannwhitneyu(plt_data[plt_data.Group==sample]['ind'], plt_data[plt_data.Group==control]['ind'], use_continuity=True, alternative='two-sided') 
    # ind_p = stats.mannwhitneyu(plt_data[plt_data.Group==sample]['ind'], plt_data[plt_data.Group==control]['ind']) 
    

    print('sample: mean ind mets', simean, '+/-', sisem, 'p-val', ind_p[1])

    # embed()
    fig, ax = plt.subplots(figsize=(20,15))
    fontsize = 30
    width = 0.35
    bar_positions = np.arange(len(diag_mets)) + .35

    ax.set_xlim(0, max(bar_positions) + 1)

    ax.bar(bar_positions, plt_data.diag, color='red', width=width, label='Diagnostic')
    ax.bar(bar_positions, plt_data.ind, color='green', bottom=plt_data.diag, width=width, label='Individual')
    # ax.axhline(y=plt_data[plt_data.Group == 'control']['total'].median(), color='b', linestyle='--', linewidth=3, label='Control Median')
    # ax.axhline(y=plt_data[plt_data.Group == 'CFS']['total'].median(), color='y', linestyle='--', linewidth=3, label='CFS Median')
    
    
    ax.set_xticks(bar_positions + width/2)
    ax.set_xticklabels(plt_data.index.tolist(), rotation=90, fontsize=fontsize*.7)
    # ax.set_yticks([0,10,20,30,40,50,60,70,80])
    # ax.set_yticklabels([0,10,20,30,40,50,60,70,80], fontsize=fontsize)
    ax.set_yticks(yaxis)
    ax.set_yticklabels(yaxis, fontsize=fontsize)
    ax.set_ylabel('Metabolites out of normal range', fontsize=fontsize)
    ax.set_title('Diag/Ind Metabolite Plot (%s n=%i; %s n=%i) Z=%0.2f, VIP=%0.2f'%(sample, sample_count, control, control_count,Z_score,VIP_score), fontsize=fontsize)
    ax.tick_params(bottom="off", top="off", left="off", right="off")
    [v.set_visible(False) for k,v in ax.spines.items()]
    
    text =  '%s: mean total mets %0.2f +/- %0.2f p-val %0.2E'%(sample,stmean, stsem, total_p[1])
    text+=  '\n%s: mean diag mets %0.2f +/- %0.2f p-val %0.2E'%(sample,sdmean, sdsem, diag_p[1])
    text+= '\n%s: mean ind mets %0.2f +/- %0.2f p-val %0.2E'%(sample,simean, sisem, ind_p[1])
    
    plt.text(1, plt_data[plt_data.Group==sample]['total'].max(),text, fontsize=fontsize*.8)
    
    text =  '%s: mean total mets %0.2f +/- %0.2f'%(control,ctmean, ctsem)
    text+=  '\n%s: mean diag mets %0.2f +/- %0.2f'%(control,cdmean, cdsem)
    text+= '\n%s: mean ind mets %0.2f +/- %0.2f'%(control,cimean, cisem)
    x = len(plt_data[plt_data.Group==sample]['diag'])+1
    y = plt_data[plt_data.Group==control]['total'].max()
    plt.text(x, y, text, fontsize=fontsize*.8)
    # embed()
    
    plt.legend(fontsize=fontsize*.8)
    # plt.show()
    plt.savefig(outfile)
    # plt.savefig(dir+'diag_met_plot_Z=%0.2f_VIP=%0.2f.svg'%(Z_score,VIP_score))
    

