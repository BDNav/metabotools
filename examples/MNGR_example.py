
from metabotools.correlation_analysis.bootstrap import run_parallel_bootstraps, analyze_bootstraps, organize_raw_data
import pandas as pd



num_bootstraps = 50
start = 6
step_size = 1

data_file = 'files/input_data.csv' # NOTE: deleted index + extra (non met) columns and saved xlsx file to csv for corr analysis
tableS1_file = 'Example_ASD_TableS1.xlsx'
db_file = r'files/MasterDatabase_v7.156.xlsx'
exp = 'ASD'
control = 'TD'


use_zscore=True
log_transform=True
method = 'pearson'


''' RUN bootstraps with replacement '''
OUTDIR='Pearson_w_replacement' # where files will be saved
run_parallel_bootstraps(OUTDIR, data_file, tableS1_file, exp, control, db_file, method, start, step_size, num_bootstraps, use_zscore=use_zscore, log_transform=log_transform, leavein=True)

''' analyze output '''
input_files = glob(f'{OUTDIR}/*_{method}.xlsx')
name = 'MNGR_with_replacement'
analyze_bootstraps(name, method, input_files, data_file, tableS1_file, db_file, transform)

collect_raw_data(OUTDIR, name, method, 'CASE')
collect_raw_data(OUTDIR, name, method, 'CONTROL')


''' RUN bootstraps without replacement '''
OUTDIR='Pearson_wo_replacement' # where files will be saved
run_parallel_bootstraps(OUTDIR, data_file, tableS1_file, exp, control, db_file, method, start, step_size, num_bootstraps, use_zscore=use_zscore, log_transform=log_transform, leavein=False)

''' analyze output '''
input_files = glob(f'{OUTDIR}/*_{method}.xlsx')
name = 'MNGR_without_replacement'
analyze_bootstraps(name, method, input_files, data_file, tableS1_file, db_file, transform)

collect_raw_data(OUTDIR, name, method, 'CASE')
collect_raw_data(OUTDIR, name, method, 'CONTROL')

