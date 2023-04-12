import pandas as pd
from IPython import embed
import pylab as plt
import seaborn as sns


import io
from contextlib import redirect_stdout
from datetime import date



def main():

    # inputfile = 'Correlation_F_suicide_plasma_CSF_10191002_DepressOnly.csv'
    # inputdir = 'test_data/'
    
    inputfile = 'Plasma_M_suicide_Sho_rap_GDF+FGF_20190926_N=93.csv'
    inputdir = 'test_data\Plasma_M/'
    outfile = 'Plasma_M_Benford.xlsx'
    
    # inputfile = 'Plasma_F_suicide_Sho_Rap_FGF+GDF_20190929_N=99.csv'
    # inputdir = 'test_data\Plasma_F/'
    # outfile = 'Plasma_F_Benford.xlsx'
    
    
    data = pd.read_csv(inputdir + inputfile, index_col=0)
    # data = pd.read_csv(r'D:\Dropbox (Personal)\naviauxlab_informatics\SyntheticDataset/Simulation Metabolomics Data_linear AUCs.csv', index_col=0) # synthetic dataset
    # calc_benford(data,1)
    calc_benford2(data,1,inputfile,outfile)
    embed()
    
# using benford_py package: https://github.com/milcent/benford_py
# good article on this: https://towardsdatascience.com/frawd-detection-using-benfords-law-python-code-9db8db474cf8

#NOTE activate py37 for this to work - graphing fxns are in the jupyter notebook...

def calc_benford2(orig_data,md_col,inputfile,outfile, name='', digs=2):
    data = orig_data[orig_data.columns[md_col:]]
    
    # get rid of any nulls:
    
    # data2=data.transpose().dropna().transpose()
    # embed()
    tmp = data.sum()
    keep = tmp[tmp>0]
    data2 = data[keep.index]
    
    total_cells = len(data)*len(data.columns)
    non_na = data.stack()
    
    # data3 = abs(data2.stack()) # non negative numbers
    data3 = non_na
    
    import benford as bf
    # f1d = bf.first_digits(data2.values, digs=2, decimals=0, show_plot=False) # digs=1 for the first digit (1-9)
    
    f = io.StringIO()
    with redirect_stdout(f):
        f1d = bf.first_digits(data3.values, digs=2, decimals=0, show_plot=False, MAD=True, KS=True, chi_square=True, confidence=95)
    out = f.getvalue()
    
    f = io.StringIO()
    with redirect_stdout(f):
        benford1dig = bf.first_digits(data3.values, digs=1, decimals=0, show_plot=False, MAD=True, KS=True, chi_square=True, confidence=95)
    out1dig = f.getvalue()
    # embed()
    
    def get_benford_stats(out):
    
        MAD = out[out.find('Mean Absolute Deviation:'):].split('\n')[0]
        MAD_c = out[out.find('MAD <='):].split('\n')[0]
        
        CS = out[out.find('The Chi-square statistic is'):].split('\n')[0]
        CS_c = out[out.find('Critical Chi-square for this series:'):].split('\n')[0]
        
        KS = out[out.find('The Kolmogorov-Smirnov statistic is '):].split('\n')[0]
        KS_c = out[out.find('Critical K-S for this series:'):].split('\n')[0]
        summary = {'summary':
        [
         MAD, MAD_c, CS, CS_c, KS, KS_c
        ]}
        return summary
        
    summary = get_benford_stats(out)
    summary1dig = get_benford_stats(out1dig)
    
    
    from scipy.stats import chisquare
    # code borrowed from: https://github.com/erdogant/benfordslaw/blob/master/benfordslaw/benfordslaw.py
    
    tstats, Praw = chisquare(f1d['Counts'], f_exp=f1d['Expected'])
    CS_p = 'The Chi-square p-value is ' + str(Praw)
    summary['summary'].append(CS_p)
    
    tstats1d, Praw1d = chisquare(benford1dig['Counts'], f_exp=benford1dig['Expected'])
    CS_p1d = 'The Chi-square p-value is ' + str(Praw1d)
    summary1dig['summary'].append(CS_p1d)
        
    summary = pd.DataFrame(summary)
    summary1dig = pd.DataFrame(summary1dig)
    
    # 2 digist
    f1d['d1'] = f1d.index.astype(str).str[0]
    f1d['d2'] = f1d.index.astype(str).str[1]
    
    #1 digits
    benford1dig['d1'] = benford1dig.index.astype(str).str[0]
    
    
    found =  pd.pivot_table(f1d, index='d1', columns='d2', values='Found')
    expected = pd.pivot_table(f1d, index='d1', columns='d2', values='Expected')
    
    dif = found/expected
    
    dif1dig =  benford1dig['Found']/ benford1dig['Expected']
    
    # sns.heatmap(dif, cmap='RdBu')
    # plt.show()
    from scipy import stats
    import numpy as np
    today = date.today()
    f1d = f1d.rename(columns={'Counts':'Observed Counts'})
    f1d['Expected Counts'] = f1d['Expected']*len(non_na)
    
    benford1dig = benford1dig.rename(columns={'Counts':'Observed Counts'})
    benford1dig['Expected Counts'] = benford1dig['Expected']*len(non_na)
    
    
    ''' OLD CODE '''
    '''
    def chi_square_test(data_count,expected_counts):
        import math
        """Return boolean on chi-square test (8 degrees of freedom & P-val=0.05)."""
        chi_square_stat = 0  # chi square test statistic
        for data, expected in zip(data_count,expected_counts):

            chi_square = math.pow(data - expected, 2)

            chi_square_stat += chi_square / expected

        print("\nChi-squared Test Statistic = {:.3f}".format(chi_square_stat))
        print("Critical value at a P-value of 0.05 is 15.51.")    
        return chi_square_stat
    chi_squared_stat = chi_square_test(f1d['Observed Counts'],f1d['Expected Counts'])

    from scipy import stats
    crit = stats.chi2.ppf(q = 0.95, # Find the critical value for 95% confidence*
                      df =len(f1d)-1) 
    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,  # Find the p-value
                             df=len(f1d)-1)
    CS_p = 'The Chi-square p-value is ' + str(p_value)
    
    
    from skgof import ks_test, cvm_test, ad_test
    # res = cvm_test(f1d['Observed Counts'].tolist(),f1d['Expected Counts'].tolist())
    
    from scipy.stats import norm, uniform
    res = cvm_test((1, 2, 3), uniform(0, 4))
    '''
    
    
    
    # embed()
    
    
    readme = {'readme': {
        'Input file':inputfile,
        'Date of Benford Analysis':today.strftime("%B %d, %Y"),
        'Number of samples (rows)':len(data),
        'Number of Cells':total_cells,
        'Number of chemicals (columns)':len(data.columns),
        'Number of N/As (NaN)':total_cells-len(non_na),
        'Number of AUCs with a number':len(data2.columns),
        # 'Number of experimental groups':len(orig_data['Group'].unique()),
        'Sample type':'mixed',
        'Minimum AUC':data3.min(),
        'Maximum AUC':data3.max(),
        'Mean':data3.mean(),
        'Median':data3.median(),
        'Geomean':stats.gmean(data3),
        # 'SD*':np.log2(abs(data3)).std(),
        'SD*':2**np.log2(data3).std(),
        'Linear Space Standard Deviation':data3.std(),
        'Interquartile range (IQR: 25th to 75th percentile)':'%i - %i '%(data3.quantile(0.25), data3.quantile(0.75)),
        '2.5th Percentile':(data3.quantile(0.025)),
        '97.5th Percentile':(data3.quantile(0.975)),
        'Skewness':data3.skew(),    
        'Kurtosis':data3.kurtosis(),
        'N=Total Numbers in Array':len(non_na),
    }}
    readme = pd.DataFrame(readme)
    
    
    
    writer = pd.ExcelWriter(outfile)
    readme.to_excel(writer,'readme')
    summary.to_excel(writer,'2d_summary')
    f1d.to_excel(writer,'2d_benford')
    dif.to_excel(writer,'2d_difference')
    
    summary1dig.to_excel(writer,'1d_summary')
    benford1dig.to_excel(writer,'1d_benford')
    dif1dig.to_excel(writer,'1d_difference')
    
    writer.save()
    
    sns.heatmap(dif, cmap='RdBu_r', vmin=0.6, vmax=1.4)
    plt.savefig('dif_plot.png')
    plt.close()
    
    f2d = bf.first_digits(data3.values, digs=2, decimals=0, show_plot=True, MAD=True, KS=True, chi_square=True, confidence=95)
    if name!='':
        plt.savefig(name+'_2D.png')
        plt.savefig(r'C:\Users\jon\pictures\%s_2D_benford_barchart_plot.png'%name)
    else:
        plt.savefig(r'C:\Users\jon\pictures\barchart_plot.png')
    f1d = bf.first_digits(data3.values, digs=1, decimals=0, show_plot=True, MAD=True, KS=True, chi_square=True, confidence=95)
    if name!='':
        plt.savefig(name+'_1D.png')
        plt.savefig(r'C:\Users\jon\pictures\%s_1D_benford_barchart_plot.png'%name)
    else:
        plt.savefig(r'C:\Users\jon\pictures\barchart_plot.png')
    plt.close()
    
'''
def calc_benford(data, md_col):
    
    data = data.astype(str)
    out = {}
    for c in data.columns:
        data[c]=data[c].str[0]
        tmp = data[c].value_counts()
        out[c]=tmp
    benford = pd.DataFrame(pd.DataFrame(out).transpose().sum())
    benford['ratio'] = benford/benford.sum()
    
    # benford.plot(kind='bar')
    # plt.show()
    
    embed()
'''



if __name__=='__main__':
    main()