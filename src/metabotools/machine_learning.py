

# comparing ML to that in metaboanalyst
from IPython import embed
import pandas as pd
import numpy as np

from correlation_analysis import correlation_analysis

dir = r'D:\Dropbox (Personal)\naviauxlab_informatics\Suicide_Study\2_26_2021\Plasma_F'
data_file = dir+'/Plasma_F_suicide_Sho_Rap_FGF+GDF_20190929_N=99.csv'
exp = 'Depress'
control = 'Control'

data = pd.read_csv(data_file,index_col=0)


zscores = correlation_analysis.get_zscores(data, exp, control)
zscores = zscores['case'].append(zscores['control'])
zscores['Group']=data['Group']
# zscores.to_csv('test_data/my_data_normalized.csv')
'''
 NOTE metaboanalyst does log10 transforms and calculates z-scores without accounting for case vs controls
'''


# norm = pd.read_csv('test_data/data_normalized.csv')
# norm.transpose().to_csv('test_data/metaboanalyst_data_normalized.csv')


# plsda on zcores
X = zscores[zscores.columns[:-1]]
y = zscores['Group']
y[y==exp]=1
y[y==control]=0
y = pd.to_numeric(y)


from sklearn.cross_decomposition import PLSRegression
plsr = PLSRegression(n_components=5, scale=False) # <1>

def calculate_vips(model):
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_
    p, h = w.shape
    vips = np.zeros((p,))
    s = np.diag(np.matmul(np.matmul(np.matmul(t.T,t),q.T), q)).reshape(h, -1)
    total_s = np.sum(s)
    for i in range(p):
        weight = np.array([ (w[i,j] / np.linalg.norm(w[:,j]))**2 for j in range(h) ])
        vips[i] = np.sqrt(p*(np.matmul(s.T, weight))/total_s)
    return vips

plsr.fit(X, y) # <2>

# vips = calculate_vips(plsr)
vips = pd.DataFrame(vips, index=X.columns)
embed()

plsr_weights = pd.DataFrame(plsr.x_weights_, index=X.columns).sort_values(0)
# embed()

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
lda = LDA(n_components=2)
Xlda = lda.fit_transform(X,y)


lda_coef = pd.DataFrame(lda.coef_, columns=X.columns).transpose()
embed()



