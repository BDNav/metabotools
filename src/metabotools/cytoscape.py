import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from IPython import embed


def get_cytoscape_output(output, db_file, VIP_score, full_data, table1_mets):
    full_data = pd.read_csv(full_data).dropna(axis=1, how='all')
    alias = pd.read_excel(db_file,'Alias')
    chem = pd.read_excel(db_file,'Chemical')
    mrm = pd.read_excel(db_file,'MRM')

    cytoscape = pd.merge(alias, table1_mets, left_on='MRM Name', right_on='MRM Name')
    cytoscape = pd.merge(cytoscape, mrm[['MRM Index', 'Chemical Index']], left_on='MRM Index', right_on='MRM Index')
    cytoscape = pd.merge(cytoscape, chem, left_on='Chemical Index', right_on='Chemical Index')

    cytoscape[['MRM Name','Z Score','Chemical Name','Cytoscape Name']].to_excel(output)

def gen_MRM_to_cytoscape_list(db_file, output_filename):
    # alias = pd.read_excel(db_file,'Alias')
    chem = pd.read_excel(db_file,'Chemical')
    mrm = pd.read_excel(db_file,'MRM')

    

    # cytoscape = pd.merge(alias, table1_mets)
    output = pd.merge(chem, mrm, left_on='Chemical Index',right_on='Chemical Index')
    # cytoscape = pd.merge(cytoscape, chem)
    output = output[['MRM Name','Chemical Name','Cytoscape Name']].dropna()
    
    output.to_excel(output_filename)
    
def find_unmapped_nodes(mapped_nodes, cytoscape_file): # compare results of cytoscape output to mapped nodes file:
    mapped_nodes = pd.read_csv(mapped_nodes,index_col=0) # nodes in the cytoscape map
    cytoscape_file = pd.read_excel(cytoscape_file) # this is made in the get_cytoscape_output() function
    # can get this file by first mapping nodes as below, then exporting table from cytoscape
    # FIRST: select all nodes in the map (use ctrl-A)
    # go to "Export Table" (figure that looks like a table going out)
    # Selected "Total Metabolism Default Node"
    # When importing be sure to set Cytoscape name as Primary Key
    # embed()
    mapped_nodes=mapped_nodes.fillna('')
    # mapped_nodes = mapped_nodes[~mapped_nodes.name.str.startswith('Node')]
    # mapped_nodes = mapped_nodes[~mapped_nodes.name.str.contains('Bond')]
    # mapped_nodes = mapped_nodes[~mapped_nodes.name.str.startswith('Saturated')]
    mapped_nodes = mapped_nodes[mapped_nodes.selected==True]
    mapped_nodes['IN CYTOSCAPE']=True
    mapped_nodes=mapped_nodes[['name','IN CYTOSCAPE']].drop_duplicates()
    print(len(mapped_nodes), 'total nodes in the map')
    
    
    out = pd.merge(mapped_nodes,cytoscape_file, left_on='name', right_on='Cytoscape Name', how='outer')
    out = out.drop_duplicates()
    # print len(out)
    out = out.fillna('-')
    # embed()
    
    
    on_map_missing = out[(out['IN CYTOSCAPE']==True)&(out['Cytoscape Name']=='-')][['name','IN CYTOSCAPE']]
    print(len(on_map_missing),)
    print('in map, missing a cytoscape mapping')
    
    
    missing_on_map = out[out['IN CYTOSCAPE']!=True][['MRM Name','Z Score','Chemical Name','Cytoscape Name']]
    print(len(missing_on_map),)
    # embed()
    print('measured, not in cytoscape map')
    
    return missing_on_map, on_map_missing
    

'''
Notes on building cytoscape map
1. export for_cytoscape figure as above
2. copy most recent cytoscape map from dB/cytoscape (currently Metabolism_WA22_ATB6_JM3_Color.cys)
3. import data from cytoscape file (if error re primary key comes up delete pandas index, resave with new name and try again)\
3.5 ** BE SURE TO SET CYTOSCAPE NAME AS PRIMARY KEY - NOT MRM NAME
4. continuous mapping to z score for each node
5. export as svg, ungroup in illustrator, select same appearance of nodes make width 0.25, 
6. add legend, 
7. export as PDF (using save as)
'''

if __name__=='__main__':
    # mapped_nodes='cytoscape_data/Metabolism_WA22_ATB6_JM5_Color_mapped_nodes.csv'
    mapped_nodes = r'D:\Dropbox (Personal)\naviauxlab_informatics\packages\metabo\cytoscape_data\Metabolism_WA22_ATB6_JM6_Color_mapped_nodes.csv'

    cytoscape_file='cytoscape_data/CFS_cytoscape.xlsx'
    out = find_unmapped_nodes(mapped_nodes, cytoscape_file)
    embed()


