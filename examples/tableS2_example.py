

import metabotools
import pandas as pd



''' SETUP '''
DIR = 'files/'
data_file = DIR+'input_data.csv'
vip_file = DIR+'plsda_vip.csv'
mda_file = DIR+'randomforests_sigfeatures.csv'
exp = 'ASD'
control = 'TD'
study_name = 'Example'


metric_cutoff=1.5

# name of the database file to use for the analysis:
db_file = 'files/MasterDatabase_v7.156.xlsx'



''' TABLE S2 '''
print('Running Table S2 correlations for %s this may take some time'%study_name)
correlation_summary_file = DIR + study_name + 'TableS2.xlsx'
from metabotools.correlation_analysis import correlation_analysis
corr_summary = correlation_analysis.correlation_summary_tables2(data_file, control, exp, db_file, correlation_summary_file) # this will now be TableS2
print('Saved TableS2 to %s'%correlation_summary_file)