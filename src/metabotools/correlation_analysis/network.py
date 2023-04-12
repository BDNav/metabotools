import pandas as pd
from IPython import embed
import pylab as plt
import numpy as np


# some nice info on adjusting the network: http://www.coppelia.io/2014/07/an-a-to-z-of-extra-features-for-the-d3-force-layout/
import os
db = ''

def main():
    # db_version = '4.4.2'
    # db_file = r'D:\Dropbox (Personal)\naviauxlab_informatics\mtdb\MasterDatabase_v%s.xlsx'%db_version
    # dir = r'D:\Dropbox (Personal)\naviauxlab_informatics\GWI/'
    # data_file= dir+'Gulf_war_reanalyzed_20170317.csv'
    # init(db_file)

    # data = pd.read_csv(data_file)
    # data.drop(data.columns[[2]], axis=1, inplace=True)
    
    # groups = ['Control','GWI']
    # group='Control'
    # method='spearman'
    # threshold=0.8

    # network = get_corr(data, group, method, threshold=threshold, Z=True, control_name='Control')
    outdir = r'D:\Dropbox (Personal)\naviauxlab_informatics\Correlation Analyses\Visualization\zscore_networks\pearson\CFS2'
    
    
    name = 'test'
    sex='male'
    sample='control'
    corr='pearson'
    in_out='all'
    generate_viz_index(outdir, sex, sample, corr, in_out, name)


    
def init(db_file):
    # dir = '\\'.join(__file__.split('\\')[:-3])
    # db_file = '/mtdb/MasterDatabase_v4.1.5.xlsx'
    global db
    
    gc = pd.read_excel(db_file, 'GC_MRM')
    lc = pd.read_excel(db_file, 'LC_MRM')
    
    shared_cols = list(set(gc.columns.tolist()) & set(lc.columns.tolist()))
    shared_cols.append('BestColumn')
    # mrm = lc[shared_cols].append(gc[shared_cols])
    mrm = lc.append(gc)
    mrm = mrm[shared_cols]
    # mrm = pd.read_excel(db_file,'MRM') # updated 12/1/2020 for new combined MRM DB
    
    # alias = pd.read_excel(db_file,'Alias', index_col=0)
    alias = pd.read_excel(db_file,'Alias')
    
    alias = alias[['MRM Name','MRM Index']].drop_duplicates()
    chem = pd.read_excel(db_file,'Chemical', index_col=0)
    pwy = pd.read_excel(db_file,'Pathway', index_col=0)
    
    # del mrm['MRM Name']
    
    db = pd.merge(mrm, alias, left_on='MRM Index', right_on='MRM Index')
    db = pd.merge(db, chem, left_on='Chemical Index', right_on='Chemical Index')
    db = pd.merge(db, pwy, left_on='Pathway Index', right_on='Pathway Index')
    
    
def get_corr(data, group, method, threshold=0.9):
    to_corr = data[data['Group']==group][data.columns[2:]]
    # corr=to_corr.corr(method=method, min_periods=1)
    import numpy as np
    
    
    corr=to_corr.corr(method=method, min_periods=1)
            
            # get all corr values as a list
    corr_list = corr.unstack()
    
    ranked_corr = corr_list.copy()
    # embed()
    ranked_corr.sort_values()
    ranked_corr=ranked_corr[ranked_corr<1]
    print('top corr', ranked_corr.idxmax())
    print('top anticorr', ranked_corr.idxmin())
            
    # plot top/bottom corrs:
    '''
    for i in ranked_corr.index[-5:]:
        
        met1 = i[0]
        met2 = i[1]
        x = to_corr[met1]
        y = to_corr[met2]
        plt.scatter(x,y)
        plt.xlabel(x_meas)
        plt.ylabel(y_meas)
        plt.show()
    '''
            
    mask = threshold
    # embed()
    masked = ranked_corr[abs(ranked_corr)>mask]
    network = []
    # build network:
    for i in masked.index:
        
        val = masked.loc[i]
        met1 = i[0]
        met2 = i[1]
        network.append({'from':met1, 'to':met2, 'corr':val})
    network = pd.DataFrame(network)
    return network

    
def analyze_islands(file_format):
    print('analyzing islands')
    print(file_format)
    from glob import glob
    files = glob('*'+file_format+'*.csv')
    all_data = pd.DataFrame()
    for f in files:
        data = pd.read_csv(f, index_col=0)
        data['source']=f.split('/')[-1].split('.')[0]
        all_data=all_data.append(data)
    all_data.to_csv('islands.csv')
    all_data['val']=1
    pivot = pd.pivot_table(all_data, values='val', index='node', columns='source')
    pivot = pivot.fillna(0)
    import seaborn as sns
    import pylab as plt
    tmp = pivot.transpose().sum()
    to_keep = tmp[tmp<4].index
    pivot = pivot.ix[to_keep]
    sns.clustermap(pivot)
    plt.show()
    embed()
    
def build_network(network, filename='d3network', type='all', set_color='group', to_color=[]):
    import networkx as nx
    G=nx.Graph()
    missing_pwys = []
    for i in network.index:
        source=network.loc[i,'from']
        target=network.loc[i,'to']
        corr = network.loc[i,'corr']
        color='lightgray'
        
        if source[-2:]=='_T':
            source = source[:-2]
        if target[-2:]=='_T':
            target = target[:-2]
        
        
        
        if corr<0:
            color='black'
            
        # error handling (for compounds with missing pathways)
        '''
        try:
            source_pwy = db[db['MRM Name']==source]['Pathway Name'].values[0]
        except:
            # embed()
            source_pwy = 'Unknown'
            missing_pwys.append(source)
        try:
            target_pwy = db[db['MRM Name']==target]['Pathway Name'].values[0]
        except:
            target_pwy = 'Unknown2'
            missing_pwys.append(target)
        '''
        
        # no error handling:
        # '''
        source_pwy = db[db['MRM Name']==source]['Pathway Name'].values[0]
        target_pwy = db[db['MRM Name']==target]['Pathway Name'].values[0]
        # '''
        
        if type=='all':
            G.add_edge(source, target, corr=corr, color=color)
            # G.add_edge(source, target, corr=corr)
        elif type=='in_pwy':
            if source_pwy==target_pwy:
                G.add_edge(source, target, corr=corr, color=color)
                # G.add_edge(source, target, corr=corr)
        elif type=='out_pwy':
            if source_pwy!=target_pwy:
                G.add_edge(source, target, corr=corr, color=color)
                # G.add_edge(source, target, corr=corr)
                
    if len(missing_pwys)>0:
        print('%i metabolites are missing pathways'%len(missing_pwys))
        print(missing_pwys)
        embed()
    for n in G.node.keys():
        G.node[n]['name']=n
    
    # get subgraphs/islands
    
    import seaborn as sns
    num_colors=20
    # palette = sns.color_palette("Set2", num_colors)
    palette = sns.color_palette("Paired", num_colors)
    # embed()
    pwy_colors = {}
    
    all_pwys = db['Pathway Name'].unique()
    all_pwys.sort()
    
    def rgb_to_hex(red, green, blue):
        """Return color as #rrggbb for the given color values."""
        return '#%02x%02x%02x' % (red, green, blue)
    x=0    
    
    for p in all_pwys:
        x = x%num_colors
        # print x
        pwy_colors[p] = rgb_to_hex(palette[x][0]*255,palette[x][1]*255,palette[x][2]*255)
        x+=1
    pwy_colors['Unknown']='lightgray'
    nodes = G.node.keys()
    # embed()
    for n in nodes:
        classification = db[db['MRM Name']==n]['Classification'].values[0]
        pathway = db[db['MRM Name']==n]['Pathway Name'].values[0]
        G.node[n]['color']=pwy_colors[pathway]
        # embed()
        if classification=='Endogenous':
            G.node[n]['type']='metabolite'
        else:
            G.node[n]['type']='toxin'
        
        
    # embed()
    ''' getting subgraph islands '''
    '''
    # embed()
    graphs = list(nx.connected_component_subgraphs(G))
    islands = []
    group=1
    for g in graphs:
        nodes = g.node.keys()
        for n in nodes:
            G.node[n]['group']=group
            # G.node[n]['color']='red'
            pwy = db[db['MRM Name_y']==n]['Pathway Name']
            if len(pwy)>0:
                pwy = pwy.values[0]
                G.node[n]['pathway']=pwy
                G.node[n]['color']=pwy_colors[pwy]
                if set_color=='select':
                    G.node[n]['color']='#d3d3d3'
                    if n in to_color:
                        G.node[n]['color']='red'
                        # embed()
            islands.append({'node':n, 'island':'island %s'%group})
        group+=1
    # embed()
    # pd.DataFrame(islands).to_csv(filename+'_islands.csv')
    '''
    # sg=G
    # embed()
    from write_json import convert_graph_for_cytoscapejs
    
    from networkx.readwrite import json_graph
    data = json_graph.node_link_data(G) # SO EASY for d3!
    
    '''
    for link in data['links']:
        source = link['source']
        embed()
        source = data['nodes'][source]['id']
        link['source'] = source
        
        target = link['target']
        target = data['nodes'][target]['id']
        link['target'] = target
        # embed()
    '''
       
    
    # data = convert_graph_for_cytoscapejs(sg)
    import json
    s = json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))
    # with open('networkx_test.json', 'w') as outfile:
    # with open('d3_test.json', 'w') as outfile:
        # outfile.write(s)
    folders = filename.split('/')[:-1]
    import os
    cdir = os.getcwd()
    for f in folders:
        cdir+='/'+f
        print(cdir)
        if not os.path.isdir(cdir):
            os.mkdir(cdir)
        
    outfile = open(filename, 'w')
    outfile.write(s)
    outfile.close()
    
    return G
    
def clean_neo4jdb(login):
    from neo4j import GraphDatabase
    
    driver = GraphDatabase.driver(login['uri'], auth=(login['user'],login['password']))
    session = driver.session()

    session.run('''MATCH (n)
    OPTIONAL MATCH (n)-[r]-()
    DELETE n,r''')
    
def load_corr_network_to_neo4j(login, project, my_network, node_info):
    from IPython import embed
    from neo4j import GraphDatabase
    
    driver = GraphDatabase.driver(login['uri'], auth=(login['user'],login['password']))
    session = driver.session()

    for i in node_info.index:
        node = node_info.loc[i,:]
        session.run("MERGE (c:Chemical {name: $name, pathway: $pathway, type: $type, group: $group})", name=node['MRM Name'], pathway=node['Pathway Name'], type=node['Exo/Endo'], group=node['Pathway Name']) # TODO add chemical name, other IDs?

    def add_edge(n1, n2, project, weight,p,fdr,q):
        session.run('''
        MATCH (c:Chemical {name: '%s'})
        MATCH (c2:Chemical {name: '%s'})
        MERGE (c)-[r:CORRELATES_WITH {weight:%0.4f, p:$p, fdr:$fdr, q:$q, project:$project}]-(c2)
        '''%(n1,n2,weight), p=p, fdr=fdr, q=q, project=project)
    
    def add_edge2(n1, n2, project, weight,p,fdr,q):
        session.run('''
        MATCH (c:Chemical {name: $n1})
        MATCH (c2:Chemical {name: $n2})
        MERGE (c)-[r:CORRELATES_WITH {weight:$weight, r:$r, p:$p, fdr:$fdr, q:$q, project:$project}]-(c2)
        
        ''', n1=n1, n2=n2, weight=abs(weight), r=weight, p=p, fdr=fdr, q=q, project=project)
    # '''    
    for i in my_network.index:
        n1 = my_network.loc[i,'from']
        n2 = my_network.loc[i,'to']
        corr = my_network.loc[i,'r']
        p = my_network.loc[i,'p']
        fdr = my_network.loc[i,'FDR']
        q = my_network.loc[i,'Q value']
        # print(p,fdr,q,p)
        # print('w: %0.4f, p:%f, fdr:%f, q:%f'%(corr,p,fdr,q))
        # embed()
        
        add_edge2(n1,n2,project,corr,p,fdr,q)
    # '''
    
def plot_corr_histograms():
    import seaborn as sns
    groups = data['Group'].unique()
    for group in groups:
        to_corr = data[data.Group==group][data.columns[2:]]
        methods = ['pearson','spearman']
        for method in methods:
            # calculat pairwise correlation
            corr=to_corr.corr(method=method, min_periods=1)
            
            # get all corr values as a list
            corr_list = corr.unstack()
            
            fig = plt.figure(1, figsize=(8, 5))
            
            corr_list.hist(bins=25)
            
            plt.title('%s %s correlation histogram'%(group, method))
            ax = fig.add_subplot(111)
            
            kurt = corr_list.kurtosis()
            avg = corr_list.mean()
            std = corr_list.std()
            skew = corr_list.skew()
            
            anno = "Mean: %0.2f\n SD: %0.2f\n Skew: %0.2f\n Kurtosis: %0.2f\n"%(avg, std, skew, kurt)
            ax.annotate(anno, xy=(min(corr_list), 12500), xycoords='data',
                xytext=(0, -50), textcoords='offset points')
            # embed()
            plt.xlabel('%s\'s r'%method.title())
            plt.ylabel('count')
            
            
            # plt.show()
            plt.savefig('%s_%s.png'%(group, method))
            plt.close()


    
def generate_viz_index(outdir, sex, sample, corr, in_out, name, forceadd_version=3.2):
    '''
    takes a directory with json networks in it and generates the index html file for viewing the networks
    
    '''
    outfile = '%s_%s_%s_%s'%(sex,sample,corr,in_out)
    # embed()
    from glob import glob
    thresholds = [0.99, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
    files = []
    for t in thresholds:
        files.append('%s_%s_%s_%0.2f_%s.json'%(sex,sample,corr,t,in_out))
    
    
    html = '''
    <!DOCTYPE html>
<!-- <script src="https://d3js.org/d3.v4.min.js"></script> -->
<script src="//d3js.org/d3.v3.min.js"></script>
<script type="text/javascript" src="https://code.jquery.com/jquery-1.7.1.min.js"></script>
<!-- <script src="jquery-3.1.1.slim.min.js"></script> -->
<!-- <script src="forceadd2.js"></script> -->
<script src="https://naviauxlab.ucsd.edu/internalapps/javascript/forceadd''' + str(forceadd_version) + '''.js"></script>
<script src="https://naviauxlab.ucsd.edu/internalapps/javascript/network_functions.js"></script>
<style>

.link {
  stroke: #ccc;
}

.node text {
  pointer-events: none;
  font: 15px sans-serif;
}

</style>
<body>
<center>

<h1>'''
    html+=name
    html+='''
</h1>
<h3><a href="#" onclick="addNetwork(network1)">Load Network</a></h3>
<input type="range" min="0.7" max="1" value="1" step="0.05" onchange="showValue(this.value)" />
<span id="range">1</span>
<br>
<svg id="graph" width="1650" height="1080"></svg>
<script>
<!-- var tmp = document.getElementById(); -->
graph = new myGraph("#graph");

var nodes = [];
var links = [];

var network1;
var network2;
var network3;
var network4;
var network5;
var network6;
var network7;


var network_dict={};

'''
    count = 1
    for t in thresholds:
        file='%s_%s_%s_%0.2f_%s.json'%(sex,sample,corr,t,in_out)
        html+="\
d3.json(\"%s\", function(json) { \
network%i = {nodes:json.nodes, links:json.links}; \
network_dict[%0.2f] = network%i; \
});\n"%(file,count,t,count)
        count+=1
    
    html+='''
var node_dict={};
var edge_dict={};

</script>


</center>
</body>
</html>
'''
    # embed()
    
    fout = open(outdir+'\\'+outfile+'.html','w')
    fout.write(html)
    fout.close()
    
    
def zscore_corrs_for_kefeng(data, case, control, type, method):
    from correlation_analysis import get_zscores, get_correlations, map_to_cytoscape, map_to_mdtb
    import datetime
    now = datetime.datetime.now()

    # if outfile_name=='':
        # outfile_name='%s_%s_%s.xlsx'%(study,type,method)
    
    SCRIPT_WD = '\\'.join(__file__.split('\\')[:-1])
    zscores = get_zscores(data, case, control)

    case_z = zscores['case']
    # data for CONTROLS:
    control_z = zscores['control']
    num_cases = len(case_z)
    num_controls = len(control_z)

    if type=='control':
        data = control_z
    if type=='case':
        data=case_z
    # embed()
    
    corrs = data.corr()
    
    return corrs

''' THIS IS DEPRECATED - use correlation_analysis.make_corr_spreadsheet instead
def make_corr_spreadsheet(data, case, control, type, method, db_file='../mtdb/MasterDatabase_v5.1.2.xlsx', infile='', outfile_name=''):
    from correlation_analysis import get_zscores, get_correlations, map_to_cytoscape, map_to_mdtb
    import datetime
    now = datetime.datetime.now()

    if outfile_name=='':
        outfile_name='%s_%s_%s.xlsx'%(study,type,method)
    
    SCRIPT_WD = '\\'.join(__file__.split('\\')[:-1])
    zscores = get_zscores(data, case, control)

    case_z = zscores['case']
    # data for CONTROLS:
    control_z = zscores['control']
    num_cases = len(case_z)
    num_controls = len(control_z)

    if type=='control':
        data = control_z
    if type=='case':
        data=case_z
    # embed()

    corrs = get_correlations(data, method=method)
    # embed()

    # db_file = '../mtdb/MasterDatabase_v5.1.2.xlsx'
    # Note this file and functions come from Suicide_Study
    cyto_xy = SCRIPT_WD+'\\Cytoscape_xy_with_ID_v2.xlsx'

    mtdb_map = map_to_mdtb(corrs, db_file)
    # embed()

    cyto_map = map_to_cytoscape(mtdb_map, cyto_xy, method)

    # embed()

    readme = {'Input File':infile,
        'Date':str(now.isoformat()),
        'Type':type,
        'Correlation Method':method,
        'Metabolites Measured':len(data.columns),
        'Total Pairwise Correlations':len(corrs),
        'Number of Cases':num_cases,
        'Number of Controls':num_controls,
        'Missing cytoscape nodes':len(cyto_map['no_map'])
        }
    readme = pd.DataFrame(pd.Series(readme), columns=[''])

    writer = pd.ExcelWriter(outfile_name)
    readme.to_excel(writer, 'README')
    data.to_excel(writer,'Z-Scores')
    mtdb_map.to_excel(writer,'Pairwise Correlations')
    cyto_map['cyto'].to_excel(writer, 'Cytoscape XY Map')
    cyto_map['no_map'].to_excel(writer, 'Missing on Cytoscape')
    writer.save()
'''
            
if __name__=='__main__':
    main()
    
    