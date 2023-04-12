
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



''' These files are generated in tableOne_example.py '''

table1_out = DIR+'%s_table_ONE_MWQ_FILTERED.xlsx'%(study_name)
table1_input = table1_out

table1_mets = pd.read_excel(table1_input,'Metabolites')
table1 = pd.read_excel(table1_input,'Table 1')

replace_chars = False # set this to True if you are using an old version of metaboanalyst

''' TABLE S1 '''

from metabotools.stats import new_build_table_s1
ts1_file = DIR + '%s_%s_TableS1.xlsx'%(study_name,exp)
print('generating Table S1: %s'%ts1_file)

new_build_table_s1(ts1_file, study_name, data_file, control, exp, table1_mets, db_file, VIP_cutoff=metric_cutoff, mda_file=mda_file, run_het=False, knn_vip_threshold=metric_cutoff, add_RF_AUCs=True, add_logistic_reg_AUCs=True, overwrite=0, replace_chars=replace_chars)

print('saving Table S1 to %s'%ts1_file)
 
