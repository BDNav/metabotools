# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Analzing metabolomics data with a variety of methods

### How do I get set up? ###

1. Clone this repo (Note all scripts tested with Python 3.8)

2. Install requirements with `pip install -r requirements.txt`

3. Install this package by navigating to the metabotools folder and install with: `pip install -e .`

4. All analyses require a MTDB database v7.156 is included in examples/files as MasterDatabase_v7.156.xlsx (others can be downloaded at http://mtdb2-env.eba-tadxr69d.us-west-2.elasticbeanstalk.com/chemicals)

### Example 1: TableOne ###

* You will need an input file in the format of examples/files/input_data.csv
* You will also need a plsda_vip.csv and randomforests_sigfeatures.csv downloaded from metaboanalyst.ca following the standard workflow
* The script examples/tableOne_example.py will run the TableOne on the input files
* Outputs: Example_table_ONE_MWQ_FILTERED.xlsx

### Example 2: TableS1 ###

* The script examples/tableS1_example.py will run the TableS1 on the input files
* This script requires the TableOne generated above as input
* Outputs: Example_ASD_TableS1.xlsx, Example_ASD_knnClusters_5.xlsx, Example_ASD_knnClusters_10.xlsx

### Example 3: TableS2 ###

* The script examples/tableS2_example.py will run the TableS2 on the input files
* Outputs: ExampleTableS2.xlsx

### Example 4: Metabolic Network Growth Rate (Vnet) Analysis� by bootstrap resamplings ###
* The script examples/MNGR_example.py will run the MNGR analysis on the input files 
* Adjust number of boostraps with num_bootstraps, starting size = 8 and step size = 2
* Outputs Folder with bootstrap data (named in file as OUTDIR), Example_Pearson_with_replacement_bootstrap_results_by_pathway_pearson.xlsx, Example_Pearson_with_replacement_bootstrap_results_pearson.xlsx, Example_with_replacement_CASE_all_sig_pairs_raw_data.xlsx, Example_with_replacement_CONTROL_all_sig_pairs_raw_data

