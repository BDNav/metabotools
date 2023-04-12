import pandas as pd
import numpy as np
# import metabotools.dbutils as dbutils
import matplotlib.pyplot as _pyplot
import scipy.stats as _stats
import numpy as _numpy

from IPython import embed
# from metabo.metabotools import db_merge_full, db_merge_mrm, db_merge_name, db_merge_path, db_merge_chem, db_read_excel, get_mrm_name



db = ''
d = ''
col = ''
g = ''
exp=''
control=''
z_cutoff=''
VIP_cutoff=''
table1 = ''
table1_mets = ''

def init(data_file, db_file, md_num, z_cutoff_l, VIP_cutoff_l, group_file, exp_l, control_l, table1_l, table1_mets_l):
    global db, col, d, g, exp, control, VIP_cutoff, z_cutoff, table1, table1_mets
    exp = exp_l
    control = control_l
    z_cutoff = z_cutoff_l
    VIP_cutoff = VIP_cutoff_l
    table1_mets = table1_mets_l
    table1 = table1_l
    
    db_dict = {'a': 'Alias',
           'c': 'Chemical',
           'm': 'MRM',
           'p': 'Pathway'}

    

    db_parts = db_read_excel(db_file, db_dict)
    db = db_merge_full(db_parts['m'], db_parts['c'], db_parts['p'])

    d = pd.read_csv(data_file).dropna(axis=1, how='all').set_index('Sample')
    # embed()
    col = d.columns.tolist()
    m = d[col[0:md_num]]

    d = np.log2(d[col[md_num:]])
    d.columns = pd.merge(pd.DataFrame(d.columns.tolist(), columns=['MRM Name']), db_parts['a'], on='MRM Name')['MRM Index']

    g = pd.read_csv(group_file).set_index('Sample') 

    db = db[db['MRM Index'].isin(d.columns.tolist())]

#########Stacked Z Score Outliers

def z_score_outliers(outfile, my_exp, my_control, show=0):
    # global db_file, db, d, col, m, d, g
    z_my_exp = pd.DataFrame(d.loc[g.index.tolist()][g.Group == my_exp].mean().sub(d.loc[g.index.tolist()][g.Group == my_control].mean()).div(d.loc[g.index.tolist()][g.Group == my_control].std()), columns=[my_exp])
    z_scores = d.loc[g.index.tolist()][g.Group == my_control].sub(d.loc[g.index.tolist()][g.Group == my_control].mean()).div(d.loc[g.index.tolist()][g.Group == my_control].std()).append(z_my_exp.T)

    # JM edit:
    z_scores = z_scores[:-1] # to remove the group name at the end of this list

    # embed()
    q_high = z_scores.apply(lambda x: len(x[x >= z_cutoff]),axis=1)
    q_low = z_scores.apply(lambda x: len(x[x <= -1*z_cutoff]),axis=1)
    # embed()
    # JMM EDIT to sort by total outliers:
    total = pd.DataFrame([q_high, q_low]).transpose()
    total['sum'] = total[0]+total[1]
    total.sort_values('sum', inplace=True, ascending=False)
    #
    # embed()

    fig = _pyplot.figure(figsize=(6.33,4.75))
    # JMM EDIT:
    # l = q_low
    # h = q_high
    # N = len(l)
    # ind = np.arange(N)

    l = total[1]
    h = total[0]
    N=len(l)
    ind = np.arange(N)

    width = 0.35
    
    # embed()

    p1 = _pyplot.bar(ind,l, color='g', width=width, label='Low') 
    p2 = _pyplot.bar(ind,h, color='r', bottom=l, width=width, label='High')

    _pyplot.ylabel('Metabolites out of normal range')
    _pyplot.xticks(ind+width/2., total.index.tolist(), rotation=90 )
    _pyplot.yticks(np.arange(0,(((q_high + q_low).max() // 25) + 3) * 25,25))
    _pyplot.xlim(0, len(ind))
    _pyplot.ylim(0, ((((q_high + q_low).max() // 25) + 2) * 25) + 5)
    _pyplot.plot([0, len(ind)], [h.loc[g[g.Group == my_control].index.tolist()].add(l.loc[g[g.Group == my_control].index.tolist()]).median(), h.loc[g[g.Group == my_control].index.tolist()].add(l.loc[g[g.Group == my_control].index.tolist()]).median()], 'y--', label='%s Median'%my_control)
    _pyplot.legend(loc=0)
    _pyplot.tight_layout()
    if show==0:
        _pyplot.savefig(outfile, dpi=600)
    else:
        _pyplot.show()

# embed()

def stacked_pathway(outfile):
    # print 'this function has been removed'
    # return
    ## REMOVING BELOW - it gives different results from HiLo plot in diagnostics, also its the same info...
    # HiLo plot gives same results from CFS1 study
    # NOTE: axes don't seem to be correcly sorted in figure
    
    z_exp = pd.DataFrame(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()), columns=[exp])
    #########Stacked horizontal pathway
    w = pd.merge(z_exp.reset_index(),db, on='MRM Index', how='left')
    e = w.groupby('Short Name')

    r = []
    for name, group in e:
        l = group[exp][group[exp] <= -1*z_cutoff].count()
        h = group[exp][group[exp] >= z_cutoff].count()
        order = group['Figure Order'].iloc[0]
        r.append([name,l,h,order])
        # embed()

    t = pd.DataFrame(r).sort_values(3, inplace=False, ascending=False)

    fig = _pyplot.figure(figsize=(8.2,6))
    t['sort'] = t[1]+t[2]
    t = t.sort_values(by='sort')

    l = t[1]
    h = t[2]
    N = len(l)
    ind = np.arange(N)
    height = 0.6
    # embed()

    p1 = _pyplot.barh(ind,l, color='g', height=height) 
    p2 = _pyplot.barh(ind,h, color='r', left=l, height=height)

    _pyplot.setp(p1, 'ec', 'None')
    _pyplot.setp(p2, 'ec', 'None')

    _pyplot.xticks(np.arange(0,55,5))
    _pyplot.yticks(ind+height/2., t[0].tolist(), fontsize=6, va='center')
    _pyplot.tick_params(
                        axis='both',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        left='off',      # ticks along the bottom edge are off
                        right='off',
                        bottom='off',      # ticks along the bottom edge are off
                        top='off',         # ticks along the top edge are off
                        )
    leg = _pyplot.legend( (p2[0], p1[0]), ('High', 'Low'), loc=0, frameon=False)
    leg = _pyplot.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
    frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

    # see text.Text, lines.Line2D, and patches.Rectangle for more info on
    # the settable properties of lines, text, and rectangles
    _pyplot.setp(ltext, fontsize='small')    # the legend text fontsize

    _pyplot.ylim(0, len(ind))
    _pyplot.xlim(0, 50)
    _pyplot.tight_layout()
    _pyplot.title('Stacked horizontal pathway Z cutoff = %0.2f'%z_cutoff)
    # _pyplot.xlabel('Abnormal Metabolites (N)')
    print( 'stacked pathway')
    # embed()
    # _pylot.show()
    _pyplot.savefig(outfile, dpi=600)
    


    
def stacked_pathway_VIP(outfile):
    #########Stacked horizontal pathway - VIP Score


    w = pd.merge(table1_mets.reset_index(),db, on='MRM Name', how='left')
    e = w.groupby('Short Name')
    # embed()

    r = []
    for name, group in e:
        # embed()
        l = group['VIP Score'][group['VIP Score'] <= -1*VIP_cutoff].count()
        h = group['VIP Score'][group['VIP Score'] >= VIP_cutoff].count()
        order = group['Figure Order'].iloc[0]
        r.append([name,l,h,order])
        # embed()

    t = pd.DataFrame(r).sort_values(3, inplace=False, ascending=False)

    fig = _pyplot.figure(figsize=(8.2,6))
    t['sort'] = t[1]+t[2]
    t = t.sort_values(by='sort')

    l = t[1]
    h = t[2]
    N = len(l)
    ind = np.arange(N)
    height = 0.6
    # embed()

    p1 = _pyplot.barh(ind,l, color='g', height=height) 
    p2 = _pyplot.barh(ind,h, color='r', left=l, height=height)

    _pyplot.setp(p1, 'ec', 'None')
    _pyplot.setp(p2, 'ec', 'None')

    _pyplot.xticks(np.arange(0,55,5))
    _pyplot.yticks(ind+height/2., t[0].tolist(), fontsize=6, va='center')
    _pyplot.tick_params(
                        axis='both',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        left='off',      # ticks along the bottom edge are off
                        right='off',
                        bottom='off',      # ticks along the bottom edge are off
                        top='off',         # ticks along the top edge are off
                        )
    leg = _pyplot.legend( (p2[0], p1[0]), ('High', 'Low'), loc=0, frameon=False)
    leg = _pyplot.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
    frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

    # see text.Text, lines.Line2D, and patches.Rectangle for more info on
    # the settable properties of lines, text, and rectangles
    _pyplot.setp(ltext, fontsize='small')    # the legend text fontsize

    _pyplot.ylim(0, len(ind))
    _pyplot.xlim(0, 50)
    _pyplot.tight_layout()
    _pyplot.title('Stacked horizontal pathway VIP cutoff = %0.2f'%VIP_cutoff)
    
    # _pyplot.xlabel('Abnormal Metabolites (N)')
    _pyplot.savefig(outfile, dpi=600)

    
##########Z Score histogram
def plot_met_z_score_histogram(outfile, lw=1., xlim=[-2,2]):
    z_exp = pd.DataFrame(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()), columns=[exp])
    data = z_exp.dropna().values.ravel()

    fig = _pyplot.figure(figsize=(7,6))
    data_mean = data.mean()
    data_std = data.std()
    high = int(data.max()) + 1 # JMM added +1 to cover full curve
    low = int(data.min()) - 2 # JMM removed -1 from inside paren and added -0.2 outside to cover full curve
    _pyplot.xlim(xlim)
    step = 0.1
    #Plot gaussian based on data mean and standard deviation
    x = np.arange(low, high + step, step)
    y = _stats.norm.pdf(x, loc = data_mean, scale = data_std)
    _pyplot.plot(x,y, color=(0.23, 0.45, 0.70, 1.00), linewidth=lw)
    
    #Plot histogram of data
    in_face = (0.70, 0.70, 0.70, 1.00)
    in_edge = (0.70, 0.70, 0.70, 1.00)
    out_color = (0.23, 0.45, 0.70)
    n, bins, patches = _pyplot.hist(data, color=out_color, bins=50, normed=True)
    _pyplot.setp(patches, 'ec', out_color)

    a = pd.DataFrame([n,bins,patches]).T
    for p in a[(a[1] <= z_cutoff) & (a[1] >= -1*z_cutoff)][2][:-1]:
        p.set_facecolor(in_face)
        p.set_edgecolor(in_edge)
    _pyplot.xlabel('Z score')
    _pyplot.ylabel('density')
    # _pyplot.title('Control' + control + 'Exp'+ exp)
    # embed()
    _pyplot.savefig(outfile, dpi=600)
    _pyplot.close()
    # _pyplot.show()
    
def jon_met_overlay_hist(outfile, lw=1, xlim=[-3,3], norm=False):
    import pylab as plt
    # g = groups
    # d = data for each met, log transformed
    
    ctrl_data = d.ix[g[g==control].dropna().index]
    
    control_mean = ctrl_data.mean()
    control_std = ctrl_data.std()
    
    exp_data = d.ix[g[g==exp].dropna().index]
    
    # z_exp = (exp_data-control_mean)/control_std
    z_exp = exp_data.sub(control_mean).div(control_std)
    # z_ctrl = (ctrl_data-control_mean)/control_std
    z_ctrl = ctrl_data.sub(control_mean).div(control_std)
    
    z_exp_unstacked = z_exp.unstack()
    z_exp_unstacked=z_exp_unstacked.clip_upper(5)
    z_exp_unstacked=z_exp_unstacked.clip_lower(-5)
    if norm==False:
        z_exp_unstacked.hist(bins=50, color='r',histtype='stepfilled', alpha=0.25)
    else:
        z_exp_unstacked.hist(bins=50, color='r',histtype='stepfilled', alpha=0.25, normed=True)
    # plt.title('exp')
    # plt.show()
    
    z_ctrl_unstacked= z_ctrl.unstack()
    z_ctrl_unstacked=z_ctrl_unstacked.clip_upper(5)
    z_ctrl_unstacked=z_ctrl_unstacked.clip_lower(-5)
    if norm==False:
        z_ctrl_unstacked.hist(bins=50, color='g',histtype='stepfilled', alpha=0.25)
    else:
        z_ctrl_unstacked.hist(bins=50, color='g',histtype='stepfilled', alpha=0.25, normed=True)
    plt.title('Green=%s (N=%i); Red=%s (N=%i)'%(control, len(g[g==control].dropna()), exp, len(g[g==exp].dropna())))
    plt.show()
    
    processed_data = pd.read_excel('z_score_histograms.xlsx','exp_z_scores')
    processed_data=processed_data.unstack()
    processed_data.hist(bins=50, color='g',histtype='stepfilled', alpha=0.25)
    plt.show()
    embed()

def plot_met_overlay_histogram(outfile, lw=1, xlim=[-3,3]):
    
    # z_exp = pd.DataFrame(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()), columns=[exp])
    # data = z_exp.dropna().values.ravel()
    
    ctrl = d.loc[g.index.tolist()][g.Group == control].sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()).stack()
    
    # JON calculate exp z scores same as control: 
    data = d.loc[g.index.tolist()][g.Group == exp].sub(d.loc[g.index.tolist()][g.Group == exp].mean()).div(d.loc[g.index.tolist()][g.Group == exp].std()).stack()
    
    
    
    data_mean = data.mean()
    data_std = data.std()
    high = int(data.max())+3 # JMM added +2 here to cover full curve...
    low = int(data.min() - 3) #JMM may need to edit this as well if curve is cutoff on left
    # embed()
    step = 0.1
    
    # embed()
    
    ctrl = ctrl.dropna()
    ctrl_mean = ctrl.mean()
    ctrl_std = ctrl.std()
    
    # embed()
    
    #Plot control
    fig = _pyplot.figure(figsize=(7,6))
    #Plot gaussian based on data mean and standard deviation
    x = np.arange(low, high + step, step)
    y = _stats.norm.pdf(x, loc = ctrl_mean, scale = ctrl_std)
    _pyplot.plot(x,y, color='g', linewidth=lw)
    
    
    #Plot histogram of data
    # FOR NORMALIZE
    n, bins, patches = _pyplot.hist(ctrl, color='g', bins=50, normed=True, histtype='stepfilled', alpha=0.25)
    
    #FOR RAW COUNTS:
    # n, bins, patches = _pyplot.hist(ctrl, color='g', bins=50, normed=False, histtype='stepfilled', alpha=0.25)
    # _pyplot.setp(patches, 'ec', 'None')
    
    
    #Plot data
    #Plot gaussian based on data mean and standard deviation
    x = np.arange(low, high + step, step)
    y = _stats.norm.pdf(x, loc = data_mean, scale = data_std)
    _pyplot.plot(x,y, color='r', linewidth=lw)
    
    #Plot histogram of data
    # NORMALIZED:
    n, bins, patches = _pyplot.hist(data, color='r', bins=50, normed=True, histtype='stepfilled', alpha=0.25)
    
    # RAW COUNTS
    # n, bins, patches = _pyplot.hist(data, color='r', bins=50, normed=False, histtype='stepfilled', alpha=0.25)
    
    _pyplot.setp(patches, 'ec', 'None')
    _pyplot.xlabel('Z score')
    # _pyplot.ylabel('density')
    _pyplot.ylabel('counts')
    _pyplot.title('Green=%s (N=%i); Red=%s (N=%i)'%(control, len(g[g==control].dropna()), exp, len(g[g==exp].dropna())))
    
    
    _pyplot.xlim(xlim)
    _pyplot.savefig(outfile, dpi=600)
    # _pyplot.show()
    
def plot_met_overlay_histogram2(outfile, lw=1, xlim=[-3,3]): # updated to show count instead of density
    # import pylab as plt
    # exp = g[g.Group==exp].index
    # ctrl = g[g.Group==control].index
    
    # embed()
    
    ctrl = d.loc[g.index.tolist()][g.Group == control].sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()).stack()
    # print len(ctrl)
    
    # ctrl.hist(bins=100, color='g')
    # plt.title('control: mean: %0.2f, std: %0.2f, count %0.2f'%(ctrl.mean(), ctrl.std(),len(ctrl)))
    # plt.show()
    
    # ctrl = d.loc[g.index.tolist()][g.Group == control]
    # ctrl = ctrl.sub(d.loc[g.index.tolist()][g.Group == control].mean())
    # ctrl = ctrl.div(d.loc[g.index.tolist()][g.Group == control].std())
    # ctrl = ctrl.stack()
    # embed()
    
    # data = d.loc[g.index.tolist()][g.Group == exp].sub(d.loc[g.index.tolist()][g.Group == exp].mean()).div(d.loc[g.index.tolist()][g.Group == exp].std()).stack()
    '''
    Please give the histograms another try use only the metabolite counts as the y-axis.  Just be sure, we are calculating the z-scores correctly, subtracting the CONTROL MEAN and deviding by the CONTROL SD.
    '''
    data = d.loc[g.index.tolist()][g.Group == exp].sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()).stack()
    # print len(data)
    # data.hist(bins=100, color='b')
    # plt.title('experiment: mean: %0.2f, std: %0.2f, count %0.2f'%(data.mean(), data.std(),len(data)))
    # plt.show()
    
    # embed()
    
    # z_exp = pd.DataFrame(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()), columns=[exp])
    # z_exp = z_exp.dropna().values.ravel()
    
    # z_control = pd.DataFrame(d.loc[g.index.tolist()][g.Group == control].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()), columns=[control])
    # z_control = z_control.dropna().values.ravel()
    
    # embed()
    
    
    # ctrl = d.loc[g.index.tolist()][g.Group == control].sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()).stack()
    
    data_mean = data.mean()
    data_std = data.std()
    high = int(ctrl.max()) # JMM added +2 here to cover full curve...
    low = int(ctrl.min()) #JMM may need to edit this as well if curve is cutoff on left
    # embed()
    step = 0.1
    
    
    ctrl = ctrl.dropna()
    ctrl_mean = ctrl.mean()
    ctrl_std = ctrl.std()
    
    
    fig, ax1 = _pyplot.subplots(figsize=(7,6))
    ax2 = ax1.twinx()
    
    
    # embed()
    #Plot histogram of data
    # n, bins, patches = _pyplot.hist(ctrl, color='g', normed=True, bins=50, histtype='stepfilled', alpha=0.25)
    n, bins, patches = ax1.hist(ctrl, color='r', bins=100, histtype='stepfilled', alpha=0.25)
    # embed()
    g_max = n.max()
    # ctrl.hist(color='g', bins=50, histtype='stepfilled', alpha=0.25)
    
    # _pyplot.setp(patches, 'ec', 'None')
    # ax2.plot(x,y, color='g', linewidth=lw)
    
    
    
    #Plot histogram of data
    n, bins, patches = ax1.hist(data, color='b', bins=100, histtype='stepfilled', alpha=0.25)
    b_max = n.max()
    ratio = b_max/g_max
    
    #Plot data
    high = int(data.max()) # JMM added +2 here to cover full curve...
    low = int(data.min())
    #Plot gaussian based on data mean and standard deviation
    x = np.arange(low, high + step, step)
    y = _stats.norm.pdf(x, loc = data_mean, scale = data_std)
    y = y*ratio
    ax2.plot(x,y, color='b', linewidth=lw)
    
    #Plot control
    high = int(ctrl.max()) # JMM added +2 here to cover full curve...
    low = int(ctrl.min())
    #Plot gaussian based on data mean and standard deviation
    x = np.arange(low, high + step, step)
    y = _stats.norm.pdf(x, loc = ctrl_mean, scale = ctrl_std)
    # y = y*ratio
    ax2.plot(x,y, color='g', linewidth=lw)
    
    left = abs(min(ctrl.min(),data.min()))
    right = max(ctrl.max(),data.max())
    
    xlim = [-1*max(left,right),max(left,right)]
    
    # n, bins, patches = _pyplot.hist(data, color='b', normed=True, bins=50, histtype='stepfilled', alpha=0.25)
    # _pyplot.setp(patches, 'ec', 'None')
    data = pd.DataFrame(data)
    # data.hist(color='b', bins=50, histtype='stepfilled', alpha=0.25)
    _pyplot.xlabel('Z score')
    ax1.set_ylabel('Count')
    ax2.set_ylabel('Density')
    _pyplot.xlim(xlim)
    
    _pyplot.title('Green=%s; Blue=%s'%(control, exp))
    # _pyplot.show()
    
    _pyplot.savefig(outfile, dpi=600)


##########Box and Whiskers
def metabolite_box_and_whiskers_plot(met_data, met, ctrl, exp_list, ctrl_agg=False, exp_agg=False, ctrl_ms=6, exp_ms=10, ax=None):
    c = ['r',
         'g',
         'y',
         'b',    
         'c',   
         'm'
         ]    
    if ax:
        _pyplot.sca(ax)
    
    grouped = met_data.groupby('Group')
    
    ctrl = grouped.get_group(ctrl)[met]
    exp = zip(exp_list, c[:len(exp_list)])
    
    # box_plot = _pyplot.boxplot(ctrl, vert=False, showfliers=False, whis=[2.5, 97.5])
    box_plot = _pyplot.boxplot(ctrl, vert=False, whis=[2.5, 97.5])
        
    if ctrl_agg:
        plot1 = _pyplot.plot(ctrl.mean(), 1, 'ko', markersize=ctrl_ms)
    else:
        plot1 = _pyplot.plot(ctrl, np.ones(len(ctrl)), 'ko', markersize=ctrl_ms)

    if exp_agg:
        plot2 = [_pyplot.plot(grouped.get_group(x[0])[met].mean(), 1, x[1] + 'o', markersize=exp_ms) for x in exp]
    else:
        plot2 = [_pyplot.plot(grouped.get_group(x[0])[met], np.ones(len(grouped.get_group(x[0])[met])), x[1] + 'o', markersize=exp_ms) for x in exp]

    _pyplot.tick_params(
                    axis='both',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='off',
                    left='off',      # ticks along the bottom edge are off
                    right='off',         # ticks along the top edge are off
                    labelleft='off') # labels along the bottom edge are off
    _pyplot.ylim(0.9, 1.1)
#    _pyplot.xlim(-abs(max(grouped[met])[1])[0] - 1, abs(max(grouped[met])[1])[0] + 1)
    _pyplot.xlabel(get_mrm_name(met, db_parts['m']), fontsize=8)

    
def all_box_whiskers(outfolder):
    z_exp = pd.DataFrame(d.loc[g.index.tolist()][g.Group == exp].mean().sub(d.loc[g.index.tolist()][g.Group == control].mean()).div(d.loc[g.index.tolist()][g.Group == control].std()), columns=[exp])
    z_exp['Abs'] = np.abs(z_exp[exp])
    z_exp.sort_values('Abs', ascending=False, inplace=True)

    # embed()
    z = z_exp[z_exp['Abs'] > z_cutoff]
    z = z.index.tolist()

    pages = int(len(z) / 20) + 1
    # embed()
    for page in range(pages):
        i = 20 * page
        u = z[i:i+20]
        fig = _pyplot.figure(figsize=(6.95,9))
        plot_data = pd.concat([pd.concat([g[g.Group == control], pd.DataFrame([exp], columns=['Group'], index=[exp])]), z_scores], axis=1)
        # embed()
        for j, x in enumerate(u):
            metabolite_box_and_whiskers_plot(plot_data, x, control, [exp], exp_agg=True, ax=fig.add_subplot(10,2,j+1))
        # _pyplot.tight_layout()
        _pyplot.savefig(dir+'{}_bw_{}.png'.format(exp,page+1), dpi=600)
        
def bubble_plot_new(outfile, table1, title, cutoff=2, adjust=True, force_text=0.1, y_cutoff=0, ymax=0, clip_upper=25, short_names=False, pwy_file='', show=False, ylog=False, xlog=False):
    path_plot = table1
    impact_col = [c for c in path_plot.columns if 'Sum' in c]
    path_plot['Impact'] = path_plot[impact_col].div(path_plot[impact_col].sum())
    path_plot['log10'] = -10 * _numpy.log10(path_plot['Hypergeometric p Value'])

    # NOTE: watch for impacts that are infinity....
    # path_plot['log10'] = path_plot['log10'].clip_upper(clip_upper) # changed to 25 on RKN request 5/31/2018
    # embed()
    path_plot['marker_size'] = path_plot['log10'].div(path_plot['log10'].max()).mul(1000) + 10
    path_plot['Impact2'] = _numpy.log10(path_plot['Impact'])
    #     path_plot = path_plot.reset_index()
    #     path_names = path_plot[(path_plot['log10'] >= 10) | (path_plot['Impact'] >= 0.05)]
    #     path_plot.to_csv('test.csv')
    fig = _pyplot.figure(figsize=(10,7.5))
    ax = fig.add_subplot(111)
    
    
    if xlog==True:
        print(path_plot['Impact'])
        path_plot['Impact'] = np.log10(1+path_plot['Impact'])
        
        # max = path_plot['Impact'].max()
        # path_plot['Impact'] = path_plot['Impact']/max
        # path_plot.loc[0,'Impact']=max
        
        print(path_plot['Impact'])
    # embed()
        
    if xlog==True:
        _pyplot.xscale("log") 
        
    _pyplot.scatter(path_plot['Impact'], path_plot['log10'], marker='o', s=path_plot['marker_size'], c=path_plot['Impact'], cmap=_pyplot.get_cmap('autumn_r'), edgecolors= "black")

    path_plot = path_plot.dropna()
    if ylog==True:
        _pyplot.yscale("log") 
    
    
    # _pyplot.ylim(bottom=0)
    # _pyplot.xlim(left=0)
    # _pyplot.show()
    
    
    
    if short_names==True: # use the "BubbleName" column in the Pathway tab of the MTDB
        pwys = pd.read_excel(pwy_file,'Pathway')
        path_plot = pd.merge(path_plot, pwys, left_on='Pathway Name', right_on='Pathway Name', how='left')
    # embed()
    tmp = path_plot.sort_values('Impact', ascending=False)
    tmp = tmp[:10]
    cutoff=0
    texts = []
    for i in tmp.index:
        x = tmp.loc[i,'Impact']
        y = tmp.loc[i,'log10']
        
        if x>cutoff:
            if short_names==True:
                text = tmp.loc[i,'BubbleName']
            else:
                text = tmp.loc[i,'Pathway Name']
        
            # ax.annotate(text, xy=(x,y), xytext=(x+0.05*x,y+0.05*y),
                    # verticalalignment='bottom')
            
            ''' USING ADJUST TEXT: '''
            # SEE HERE: https://github.com/Phlya/adjustText/blob/master/examples/Examples.ipynb
            # texts.append(ax.annotate(text, xy=(x,y), xytext=(x+0.05*x,y+0.05*y), horizontalalignment='right'))
            texts.append(ax.text(x+.01, y, text))
                   
    # _pyplot.show()
    #     for i in range(len(path_names)):
    #         _pyplot.text(path_names['Impact'].iloc(i), path_names['log10'].iloc(i), path_names['Pathway Name'].iloc(i))
    #_pyplot.xscale('log')
    # _pyplot.title(str(adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))))
    
    ''' ADJUST TEXT:'''
    
    # adjust_text(texts, save_steps=True, force_text=0.05)
    if adjust==True:
        from adjustText import adjust_text
        # adjust_text(texts, save_steps=False, force_text=force_text)
        adjust_text(texts, save_steps=True, force_text=force_text, save_format='png')
    
    
    _pyplot.gca().set_xlim(left=0)    
    _pyplot.gca().set_ylim(bottom=5) # CHANGED to 5 on RKN request 5/31/2018)
    # _pyplot.xlim([0,max(path_plot['Impact'])+0.1]) 
    _pyplot.ylabel('-10*log(pvalue)')
    _pyplot.xlabel('Fractional Impact') 
    # embed()
    # title = '%s (N=%i); %s (N=%i)'%(control, len(g[g==control]['Group'].dropna()), exp, len(g[g==exp]['Group'].dropna()))
    _pyplot.title(title)
    if ymax==0:
        ymax = path_plot['log10'].max()+1
    _pyplot.ylim([y_cutoff,ymax])
    # embed()
    # try:
    
    
    if show==False:
        _pyplot.savefig(outfile+'.png')
        _pyplot.savefig(outfile+'.svg')
    else:
    # except:
        # print 'strange error with saving these files'
        _pyplot.show()

def bubble_plot(outfile, cutoff=2, adjust=True, force_text=0.1, y_cutoff=0, ymax=0, clip_upper=25, short_names=False, pwy_file='', show=False):
    path_plot = table1
    impact_col = [c for c in path_plot.columns if 'Sum' in c]
    path_plot['Impact'] = path_plot[impact_col].div(path_plot[impact_col].sum())
    path_plot['log10'] = -10 * _numpy.log10(path_plot['Hypergeometric p Value'])

    # NOTE: watch for impacts that are infinity....
    # path_plot['log10'] = path_plot['log10'].clip_upper(clip_upper) # changed to 25 on RKN request 5/31/2018
    # embed()
    path_plot['marker_size'] = path_plot['log10'].div(path_plot['log10'].max()).mul(1000) + 10
    path_plot['Impact2'] = _numpy.log10(path_plot['Impact'])
    #     path_plot = path_plot.reset_index()
    #     path_names = path_plot[(path_plot['log10'] >= 10) | (path_plot['Impact'] >= 0.05)]
    #     path_plot.to_csv('test.csv')
    fig = _pyplot.figure(figsize=(10,7.5))
    ax = fig.add_subplot(111)
    
    _pyplot.scatter(path_plot['Impact'], path_plot['log10'], marker='o', s=path_plot['marker_size'], c=path_plot['Impact'], cmap=_pyplot.get_cmap('autumn_r'), edgecolors= "black")

    path_plot = path_plot.dropna()
    
    from adjustText import adjust_text
    
    if short_names==True: # use the "BubbleName" column in the Pathway tab of the MTDB
        pwys = pd.read_excel(pwy_file,'Pathway')
        path_plot = pd.merge(path_plot, pwys, left_on='Pathway Name', right_on='Pathway Name', how='left')
    # embed()
    texts = []
    for i in path_plot.index:
        x = path_plot.loc[i,'Impact']
        y = path_plot.loc[i,'log10']
        
        if x*y>cutoff:
            if short_names==True:
                text = path_plot.loc[i,'BubbleName']
            else:
                text = path_plot.loc[i,'Pathway Name']
        
            # ax.annotate(text, xy=(x,y), xytext=(x+0.05*x,y+0.05*y),
                    # verticalalignment='bottom')
            
            ''' USING ADJUST TEXT: '''
            # SEE HERE: https://github.com/Phlya/adjustText/blob/master/examples/Examples.ipynb
            # texts.append(ax.annotate(text, xy=(x,y), xytext=(x+0.05*x,y+0.05*y), horizontalalignment='right'))
            texts.append(ax.text(x+.01, y, text))
                   
    # _pyplot.show()
    #     for i in range(len(path_names)):
    #         _pyplot.text(path_names['Impact'].iloc(i), path_names['log10'].iloc(i), path_names['Pathway Name'].iloc(i))
    #_pyplot.xscale('log')
    # _pyplot.title(str(adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))))
    
    ''' ADJUST TEXT:'''
    # adjust_text(texts, save_steps=True, force_text=0.05)
    if adjust==True:
        # adjust_text(texts, save_steps=False, force_text=force_text)
        adjust_text(texts, save_steps=True, force_text=force_text, save_format='png')
    
    
    _pyplot.gca().set_xlim(left=0)    
    _pyplot.gca().set_ylim(bottom=5) # CHANGED to 5 on RKN request 5/31/2018)
    # _pyplot.xlim([0,max(path_plot['Impact'])+0.1]) 
    _pyplot.ylabel('-10*log(pvalue)')
    _pyplot.xlabel('Fractional Impact') 
    # embed()
    _pyplot.title('%s (N=%i); %s (N=%i)'%(control, len(g[g==control]['Group'].dropna()), exp, len(g[g==exp]['Group'].dropna())))
    if ymax==0:
        ymax = path_plot['log10'].max()+1
    _pyplot.ylim([y_cutoff,ymax])
    # embed()
    # try:
    if show==False:
        _pyplot.savefig(outfile+'.png')
        _pyplot.savefig(outfile+'.svg')
    else:
    # except:
        # print 'strange error with saving these files'
        _pyplot.show()
    
def bubble_plot_web(table1, cutoff=2, adjust=True, force_text=0.1, y_cutoff=0, ymax=10):
    path_plot = table1
    # embed()
    path_plot['Impact'] = path_plot['Impact (Sum Z Score)'].div(path_plot['Impact (Sum Z Score)'].sum())
    path_plot['log10'] = -10 * _numpy.log10(path_plot['Hypergeometric p Value'])

    path_plot['log10'] = path_plot['log10'].clip_upper(25) # changed to 25 on RKN request 5/31/2018
    
    path_plot['marker_size'] = path_plot['log10'].div(path_plot['log10'].max()).mul(1000) + 10
    path_plot['Impact2'] = _numpy.log10(path_plot['Impact'])
       
    out = path_plot[['Impact','log10','marker_size','Pathway Name']]
    out = out.dropna()
    out = out.rename(columns={'Impact':'x','log10':'y','marker_size':'v','Pathway Name':'name'})
    # _pyplot.scatter(path_plot['Impact'], path_plot['log10'], marker='o', s=path_plot['marker_size'], c=path_plot['Impact'], cmap=_pyplot.get_cmap('autumn_r'))
    out_json = []
    for i in out.index:
        x = out.loc[i,'x']
        y = out.loc[i,'y']
        v = out.loc[i,'v']
        name = out.loc[i,'name']
        out_json.append({'x':x,'y':y,'v':v,'name':name})

    return out_json
