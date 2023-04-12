

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


at = 'v' # the analysis type, v=VIP, z=Z score

# VIP or Z score cutoff:
z_cutoff = 2

# Q value cutoff
Q_cutoff = 1 # set to 1 to bypass this cutoff

# U value cutoff
U_cutoff = 1 # set to 1 to bypass this cutoff

# MDA cutoff
mda_cutoff = None # set to None to bypass this cutoff

replace_chars = False # set this to True if you are using an old version of metaboanalyst


""" TABLE ONE """
from metabotools.tableOne import buildTableOne
md_num=1 # number of metadata columns
table1_out = DIR+'%s_table_ONE_MWQ_FILTERED.xlsx'%(study_name)
table1_input = table1_out

print('Building Table 1 ...')
buildTableOne(data_file, exp, control, vip_file, md_num, at, metric_cutoff, db_file, table1_out, Q_filter_val=Q_cutoff, U_filter_val=U_cutoff, mda_filter=mda_cutoff, mda_file=mda_file, replace_chars=replace_chars)
print('Saving Table1 to: ', table1_out)