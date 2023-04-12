import os
import pandas as pd
from IPython import embed
import seaborn as sns

#CIRCOS_PATH = r'C:\Users\Jon\Documents\BITBUCKET\PyCircos\circos-0.69-7\bin\circos' # jons old pc
CIRCOS_PATH = r'C:\circos-0.69-9\bin\circos'

cmd ='perl %s'%CIRCOS_PATH

# os.system(cmd + ' -modules') # show modules



def main():
    mtdb = '../../../../mtdb/MTDB_v7.391.xlsx'
    
    case = pd.read_excel('../Plasma_F_Suicide_April2022TableS2.xlsx', 'pearson_case', index_col=0)
    gen_circos(case, 'F_Suicide_circos_case_July20_2', mtdb)
    
    control = pd.read_excel('../Plasma_F_Suicide_April2022TableS2.xlsx', 'pearson_control', index_col=0)
    gen_circos(control, 'F_Suicide_circos_control_July20_2', mtdb)
    
def gen_circos(data, name, mtdb, qfilter=1, rfilter=0, DIR='', posneg=None):
    copy_circos_conf()
    
    if posneg=='POS':
        data = data[data['r']>0]
    if posneg=='NEG':
        data = data[data['r']<0]
    
    import numpy as np
    
    met2path = data[['from','Pathway1']].drop_duplicates()
    met2path = met2path.rename(columns={'from':'met','Pathway1':'pathway'})
    met2path2 = data[['to','Pathway2']].drop_duplicates()
    met2path2 = met2path2.rename(columns={'to':'met','Pathway2':'pathway'})
    met2path = pd.concat([met2path, met2path2]).drop_duplicates()

    met2path = met2path.sort_values(['pathway','met'])
    
    '''
    order = met2path.groupby('pathway').count().sort_values('met')
    half = int(len(order)/2)
    bottom = order[:half] # smallest pwys
    top = order[half:] # largest pwys
    top = top.sort_values('met',ascending=False)
    new_order = []
    for i in range(0,len(top),1):
        if i < len(top):
            new_order.append(top.index[i])
        if i < len(bottom):
            new_order.append(bottom.index[i])
    '''
    
    ''' move ceramides '''
    # new_order.pop(new_order.index('Ceramide Metabolism'))
    # new_order.insert(40,'Ceramide Metabolism')
    
    # new_order.pop(new_order.index('Deoxysphingolipids'))
    # new_order.insert(39,'Deoxysphingolipids')
    
    # new_order.pop(new_order.index('Glycosphingolipid Metabolism'))
    # new_order.insert(38,'Glycosphingolipid Metabolism')
    # if 'Protein' in new_order:
    #     new_order.pop(new_order.index('Protein')) # removing this here removes it from the circos diagram
    
    
    # TODO: consider specifying full order here - remove order, half, top, etc above
       

    ''' Only plot Out of pathway correlations '''
    data = data[data['Pathway']=='OUT']
    data = data[data['Q value']<qfilter]
    data = data[data['abs(r)']>rfilter]
        
    pos_cnt = len(data[data['r']>0])
    neg_cnt = len(data[data['r']<0])
    
    filtered_corrs = len(data)
    print(filtered_corrs, ' out of pathway correlations')

    mapping = []

    
    
    pwy_data = pd.read_excel(mtdb, 'Pathway')
    pwy_data[['Pathway Name','Short Name', 'BubbleName']]
    pwy_data.index = pwy_data['Pathway Name']
    short_pwy_names = pwy_data['Short Name'].to_dict()
    
    
    ''' updating with bob's latest short names '''
    for pwy_name in PWY_NAME_UPDATES.keys():
        short_name = PWY_NAME_UPDATES[pwy_name]
        if pwy_name in short_pwy_names.keys():
            short_pwy_names[pwy_name] = short_name

    # update_names = pd.read_csv('pathway_short_names.csv')
    # for i in update_names.index:
    #     pwy_name = update_names.loc[i,'name']
    #     short_name = update_names.loc[i,'short name2']
    #     if pwy_name in short_pwy_names.keys():
    #         short_pwy_names[pwy_name] = short_name
            # print('updating %s to %s'%(pwy_name, short_name))
    # embed()
    
    
    ''' modify short pwy names: '''
    short_pwy_names['Vitamin B7 (Biotin) Metabolism'] = "Biotin"
    short_pwy_names['Vitamin B5 (Pantothenate, CoA) Metabolism'] = "Pantothenate"
    
    # merge pathways
    to_merge = {}
    to_merge['Antibiotics, Antifungals, and Antiparasitics'] = "Anthropogenic chemicals"
    to_merge['Colors and dyes'] = "Anthropogenic chemicals"
    to_merge['Flame retardants'] = "Anthropogenic chemicals"
    to_merge['Food Sources, Additives, Preservatives, Colorings, and Dyes'] = "Anthropogenic chemicals"
    to_merge['Food additives and preservatives'] = "Anthropogenic chemicals"
    to_merge['Fungicides'] = "Anthropogenic chemicals"
    to_merge['Herbicides'] = "Anthropogenic chemicals"
    to_merge['Industrial byproducts'] = "Anthropogenic chemicals"
    to_merge['OTC and Prescription Pharmaceutical Metabolism'] = "Anthropogenic chemicals"
    to_merge['Persistent organic pollutant (POPs)'] = "Anthropogenic chemicals"
    to_merge['Personal care products'] = "Anthropogenic chemicals"
    to_merge['Pesticides'] = "Anthropogenic chemicals"
    to_merge['Pesticides and Xenobiotic'] = "Anthropogenic chemicals"
    to_merge['Pharmaceuticals and veterinary drugs'] = "Anthropogenic chemicals"
    to_merge['Phthalates and plasticizers'] = "Anthropogenic chemicals"
    to_merge['Plant growth regulators'] = "Anthropogenic chemicals"
    to_merge['Surfactants'] = "Anthropogenic chemicals"
    for path in to_merge.keys():
        mets = met2path[met2path['pathway']==path]
        print(path, len(mets))
        met2path.loc[mets.index,'pathway'] = to_merge[path]

    
    # creating even shorter names:
    short_pwy_names[''] = ""
    
    global PATHWAY_ORDER

    pathways = met2path['pathway'].unique()
    colors = sns.color_palette("hls", len(pathways))

    missing_pathway_order = set(pathways)-set(PATHWAY_ORDER)
    if len(missing_pathway_order)>0:
        print('error - need to add some pwys to order:')
        print(missing_pathway_order)
        embed()

    pwys = pd.DataFrame(pathways, columns=['pwy'])
    pwys.index=pwys['pwy']

    pwy_order = pd.DataFrame(PATHWAY_ORDER, columns=['pwy'])
    pwy_order.index=pwy_order['pwy']

    FINAL_PWY_ORDER = pd.merge(pwy_order, pwys, left_index=True, right_index=True).index.tolist() # w

    rgb_colors = [','.join([str(int(c[0]*255)),str(int(c[1]*255)),str(int(c[2]*255))]) for c in colors]
    # colormap = dict(zip(pathways, rgb_colors))
    
    # PATHWAY_ORDER = list(set(pathways)&set(PATHWAY_ORDER)) # only use pathways we've measured
    colormap = dict(zip(pathways, rgb_colors))
    
    ''' GEN KARYOTYPE '''

    fout=open('karyotype.txt','w')
    fout2=open('pwy_order.tsv','w')
    
    id = 1
    for pathway in FINAL_PWY_ORDER:
            
        mets = met2path[met2path['pathway']==pathway]
        # embed()
        color = colormap[pathway]
        
        count = 0
        for met in mets['met'].values:
            mapping.append({'pathway':id,'met':met,'num':count, 'pwy_name':pathway})
            count+=1
        
        start = 0
        end = len(mets)
        if pathway in short_pwy_names.keys():
            short_name = short_pwy_names[pathway]
            # pathway = pathway.replace(' ','_')
            # fout.write('chr - %i %s %i %i %s\n'%(id, short_name, start, end, color))
        else:
            print('pwy missing',pathway)
            short_name = pathway
        if short_name is np.nan:
            short_name = pathway
        short_name = short_name.replace(' ','_')
        
        fout.write('chr - %i %s %i %i %s\n'%(id, short_name, start, end, color))
        for met in mets['met'].values:
            fout2.write(f'{id}\t{short_name}\t{met}\n')
        
        id+=1
    fout.close()
    fout2.close()
    # embed()
    mapping = pd.DataFrame(mapping)

    data = pd.merge(data, mapping, left_on='from', right_on='met')
    data2 = pd.merge(data, mapping, left_on='to', right_on='met')

    links = data2[['pathway_x','num_x','pathway_y','num_y','r']]
    # embed()

    fout=open('links2.txt','w')
    for i in links.index:
        chr1 = links.loc[i,'pathway_x']
        chr2 = links.loc[i,'pathway_y']
        met1 = links.loc[i,'num_x']
        met2 = links.loc[i,'num_y']
        # fout.write('1 100 200 2 250 300 color=blue')
        r = links.loc[i,'r']
        if r > 0:
            fout.write('%i %i %i %i %i %i color=red\n'%(chr1, met1, met1, chr2, met2, met2))
        else:
            fout.write('%i %i %i %i %i %i color=vdblue\n'%(chr1, met1, met1, chr2, met2, met2))
    fout.close()
    


    ''' GEN CONF '''


    global CIRCOS_PATH

    cmd ='perl %s'%CIRCOS_PATH
    os.system(cmd)

    name = DIR + name + '_r=%0.4f_q=%0.7f'%(rfilter,qfilter) + '_%i_OOP_corrs'%filtered_corrs
    if posneg=='POS':
        name+='_POS'
    if posneg=='NEG':
        name+='_NEG'
    
    cmd = 'copy circos.png "%s.png"'%name
    os.system(cmd)
    print('copied to %s'%name)
    cmd = 'copy circos.svg "%s.svg"'%name
    os.system(cmd)
    print('copied to %s'%name)

    return filtered_corrs

if __name__== '__main__':
    main()

global PWY_NAME_UPDATES, PATHWAY_ORDER

PWY_NAME_UPDATES = {
    '1-Carbon, Folate Metabolism': '1-Carbon/Folate',
    'Amino Acid Metabolism (not otherwise covered)': 'Thr/Asn',
    'Amino-Sugar, Galactose, & Non-Glucose Metabolism': 'Nonglucose sugars',
    'Antibiotics, Antifungals, and Antiparasitics': 'Antimicrobials',
    'Bile Salt Metabolism': 'Bile Salts',
    'Bioamines and Neurotransmitter Metabolism': 'Neurotransmitters',
    'Biopterin, Neopterin, Molybdopterin Metabolism': 'Biopterin (BH4)',
    'Branch Chain Amino Acid Metabolism': 'Ile/Val/Thr/Met',
    'Cardiolipin Metabolism': 'Cardiolipin',
    'Ceramide Metabolism': 'Ceramides',
    'Cholesterol, Cortisol, Non-Gonadal Steroid Metabolism': 'Cholesterol',
    'Collagen Metabolism': 'Collagen',
    'Colors and dyes': 'Foods and additives',
    'Deoxysphingolipids': 'Deoxysphingolipids',
    'Drugs of Abuse': 'Drugs of Abuse',
    'Eicosanoid and Resolvin Metabolism': 'Eicosanoids',
    'Endocannabinoid Metabolism': 'Endocannabinoids',
    'Fatty Acid Oxidation and Synthesis': 'Acylcarnitines',
    'Flame retardants': 'Flame retardants',
    'Food Sources, Additives, Preservatives, Colorings, and Dyes': 'Foods and additives',
    'Food additives and preservatives': 'Food additives and preservatives',
    'Fungicides': 'Fungicides',
    'GABA, Glutamate, Arginine, Ornithine, Proline Metabolism': 'GABA, gluatamate',
    'Gamma-Glutamyl and other Dipeptides': 'Dipeptides',
    'Ganglioside Metabolism': 'Gangliosides',
    'Glucuronic Acid Pathway': 'Glucuronic Acid',
    'Glycolysis and Gluconeogenesis Metabolism': 'Glycolysis',
    'Glycosphingolipid Metabolism': 'Glycosphingolipids',
    'Gonadal Steroids': 'Gonadal Steroids',
    'Heme and Porphyrin Metabolism': 'Porphyrins',
    'Herbicides': 'Herbicides',
    'Herbicides ': 'Herbicides ',
    'Histidine, Histamine, Carnosine Metabolism': 'Histidine/Histamine',
    'Industrial byproducts': 'Industrial byproducts',
    'Isoleucine, Valine, Threonine, or Methionine Metabolism': 'Ile/Val/Thr/Met',
    'Ketone Body Metabolism': 'Ketones',
    'Krebs Cycle': 'Krebs Cycle',
    'Lysine Metabolism': 'Lysine',
    'MISSING': 'MISSING',
    'Microbial metabolites': 'Microbiome',
    'Microbial metabolites ': 'Microbiome',
    'Microbiome Metabolism': 'Microbiome',
    'Neuro- and Gastropeptide Hormones': 'Neuro- and Gastropeptide Hormones',
    'Nitric Oxide, Superoxide, Peroxide Metabolism': 'NO, Lipoic acid',
    'OTC and Prescription Pharmaceutical Metabolism': 'Pharmaceuticals',
    'Oligoadenylates (2-5A), Specialized Nucleotide Metabolism': 'Oligoadenylates (2-5A), Specialized Nucleotide',
    'Pentose Phosphate, Gluconate Metabolism': 'PPP, Gluconate',
    'Peroxisomal Metabolism': 'Peroxisomal Metabolism',
    'Persistent organic pollutant (POPs)': 'Persistent organic pollutant (POPs)',
    'Personal care products': 'Personal care products',
    'Pesticides': 'Pesticides',
    'Pesticides and Xenobiotic Metabolism': 'Pesticides and Xenobiotic',
    'Pharmaceuticals and veterinary drugs': 'Pharmaceuticals',
    'Phosphate and Pyrophosphate Metabolism': 'Pyrophosphate',
    'Phospholipid (PA) Metabolism': 'Phospholipid (PA)',
    'Phospholipid (PC) Metabolism': 'Phospholipid (PC)',
    'Phospholipid (PE) Metabolism': 'Phospholipid (PE)',
    'Phospholipid (PG) Metabolism': 'Phospholipid (PG)',
    'Phospholipid (PI) Metabolism': 'Phospholipid (PI)',
    'Phospholipid (PS) Metabolism': 'Phospholipid (PS)',
    'Phospholipid Metabolism': 'Phospholipids',
    'Phthalates and plasticizers': 'Plasticizers',
    'Phytanic, Branch, Odd Chain Fatty Acid Metabolism': 'Phytanic, Branch, Odd Chain Fatty Acid',
    'Phytonutrients, Bioactive Botanical Metabolites': 'Phytonutrients',
    'Plant growth regulators': 'Plant growth regulators',
    'Plasmalogen Metabolism': 'Plasmalogens',
    'Plastics, Phthalates, Parabens, and Personal Care Products': 'Plastics, Phthalates, Parabens, and Personal Care Products',
    'Polyamine Metabolism': 'Polyamine',
    'Purine Metabolism': 'Purines',
    'Pyrimidine Metabolism': 'Pyrimidines',
    'SAM, SAH, Methionine, Cysteine, Glutathione Metabolism': 'Met/Cys/GSH',
    'Sphingolipid Metabolism': 'Sphingolipids',
    'Sphingomyelin Metabolism': 'Sphingomyelins',
    'Stable Isotope Labeled Internal Standards': 'Miscelaneous',
    'Surfactants': 'Surfactants',
    'Taurine, Hypotaurine Metabolism': 'Taurine',
    'Thyroid Hormone Metabolism': 'Thyroid ',
    'Triacylglycerol Metabolism': 'Triacylglycerol',
    'Tryptophan, Kynurenine, Serotonin, Melatonin Metabolism': 'Trp/5HT',
    'Tyrosine and Phenylalanine Metabolism': 'Tyr/Phe',
    'Ubiquinone and Dolichol Metabolism': 'CoQ10',
    'Urea Cycle': 'Urea Cycle',
    'Very Long Chain Fatty Acid Oxidation': 'C22-C28',
    'Vitamin A (Retinol), Carotenoid Metabolism': 'Vitamin A',
    'Vitamin B1 (Thiamine) Metabolism': 'Vitamin B1',
    'Vitamin B12 (Cobalamin)  Metabolism': 'Vitamin B12',
    'Vitamin B2 (Riboflavin) Metabolism': 'Vitamin B2',
    'Vitamin B3 (Niacin, NAD+) Metabolism': 'Niacin',
    'Vitamin B5 (Pantothenate, CoA) Metabolism': 'Pantothenate',
    'Vitamin B6 (Pyridoxine) Metabolism': 'Vitamin B6',
    'Vitamin B7 (Biotin) Metabolism': 'Biotin',
    'Vitamin C (Ascorbate) Metabolism': 'Vitamin C',
    'Vitamin D (Calciferol) Metabolism': 'Vitamin D',
    'Vitamin E (Tocopherol) Metabolism': 'Vitamin E',
    'Vitamin K (Menaquinone) Metabolism': 'Vitamin K'
 }


PATHWAY_ORDER = [
    u'Anthropogenic chemicals',
    u'Vitamin E (Tocopherol) Metabolism',
    u'Sphingomyelin Metabolism',
    u'Phosphate and Pyrophosphate Metabolism',
    u'Ketone Body Metabolism',
    u'Fatty Acid Oxidation and Synthesis',
    u'Thyroid Hormone Metabolism',
    u'Gonadal Steroids',
    u'Purine Metabolism',
    u'Gamma-Glutamyl and other Dipeptides',
    u'Microbiome Metabolism',
    u'Protein',
    u'Phospholipid Metabolism',
    u'Vitamin B12 (Cobalamin)  Metabolism',
    u'Glucuronic Acid Pathway',
    u'Eicosanoid and Resolvin Metabolism',
    u'Vitamin B5 (Pantothenate, CoA) Metabolism',
    u'Vitamin A (Retinol), Carotenoid Metabolism',
    u'Pyrimidine Metabolism',
    u'Vitamin B7 (Biotin) Metabolism',
    u'Cholesterol, Cortisol, Non-Gonadal Steroid Metabolism',
    u'Vitamin C (Ascorbate) Metabolism',
    u'Glycolysis and Gluconeogenesis Metabolism',
    u'Biopterin, Neopterin, Molybdopterin Metabolism',
    u'Taurine, Hypotaurine Metabolism',
    u'SAM, SAH, Methionine, Cysteine, Glutathione Metabolism',
    u'Vitamin B6 (Pyridoxine) Metabolism',
    u'Branch Chain Amino Acid Metabolism',
    u'Ubiquinone and Dolichol Metabolism',
    u'Krebs Cycle',
    u'Vitamin D (Calciferol) Metabolism',
    u'Bioamines and Neurotransmitter Metabolism',
    u'Vitamin B2 (Riboflavin) Metabolism',
    u'Cardiolipin Metabolism',
    u'Peroxisomal Metabolism',
    u'Plasmalogen Metabolism',
    u'Amino Acid Metabolism (not otherwise covered)',
    u'Bile Salt Metabolism',
    u'Phytonutrients, Bioactive Botanical Metabolites',
    u'GABA, Glutamate, Arginine, Ornithine, Proline Metabolism',
    u'Endocannabinoid Metabolism',
    u'Polyamine Metabolism',
    u'Glycosphingolipid Metabolism',
    u'Deoxysphingolipids',
    u'Ceramide Metabolism',
    u'Vitamin B1 (Thiamine) Metabolism',
    u'Tyrosine and Phenylalanine Metabolism',
    u'1-Carbon, Folate Metabolism',
    u'Collagen Metabolism',
    u'Pentose Phosphate, Gluconate Metabolism',
    u'Amino-Sugar, Galactose, & Non-Glucose Metabolism',
    u'Nitric Oxide, Superoxide, Peroxide Metabolism',
    u'Tryptophan, Kynurenine, Serotonin, Melatonin Metabolism',
    u'Histidine, Histamine, Carnosine Metabolism',
    u'Urea Cycle',
    u'Vitamin B3 (Niacin, NAD+) Metabolism',
    u'Ganglioside Metabolism',
    u'Heme and Porphyrin Metabolism'
]


def copy_circos_conf():
    '''
    Writes the circos conf file
    This defines how the map will be drawn
    It must be in the directory where circos will be called
    '''
    import os
    circos_conf_loc = (os.sep).join(os.path.abspath(__file__).split(os.sep)[:-1])+os.sep+'circos.conf'
    cmd = 'copy "%s" circos.conf'%circos_conf_loc
    os.system(cmd)