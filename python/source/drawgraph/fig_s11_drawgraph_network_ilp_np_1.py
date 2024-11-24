import loompy
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import desc
import scipy
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns

import matplotlib.font_manager as fm
import networkx as nx

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import matplotlib.colors as colors
import matplotlib.cm as cmx

from matplotlib.pyplot import rc_context

#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns


#tf와 np 클러스터

filepath="D:/project/python/bioinfo/data/loom"
#filepath="D:/project/python/bioinfo/data/h5ad"
filepath="D:/project/python/bioinfo/data/graph"
filepath="D:/project/python/bioinfo/data/graph/tf/ilp"
hg19="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/hg19/"
chgnew="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/chgnew"
pheromone="D:/project/python/bioinfo/data/graph/creatematrix/pheromone/"
creatematrixpath="D:/project/python/bioinfo/data/graph/creatematrix"
chemopath="D:/project/python/bioinfo/data/graph/chemo"
tfpath="D:/project/python/bioinfo/data/graph/tf"

filename="r_fca_biohub_testis_10x.loom"
filename="s_fca_biohub_antenna_10x.loom"
filename="s_fca_biohub_head_10x.loom"
filename="Chemosensory.txt"
filename="np_10_regulon_tissue_ann_tfgene_reggene_list.txt"
filenamenp="np_ilp.txt"
filenamenpr="npr_ilp.txt"
gfilename="genes.tsv"
bfilename="barcodes.tsv"
mfilename="matrix.mtx"
annofilename="pheromone_ann.tsv"
#filename="s_fca_biohub_head_10x.h5ad"
#filename="s_fca_biohub_head_10x.loom"
#filename="s_fca_biohub_head_10x_.loom"
#filename="fca_biohub_fat_body_ss2.loom"
#filename="s_fca_biohub_head_10x.h5ad"
#filename="10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5"
results_dir = 'D:/project/python/bioinfo/data/h5ad'  # the file that will store the analysis results
results_file="graph_5gene.h5ad"


filedir= os.path.join(filepath,filenamenp)

print(filedir)

#pbarcodes = pd.DataFrame(os.path.join(chemopath,filename),index=0)
dtfm = pd.read_csv(os.path.join(filepath,filenamenp), delimiter='\t')
#genename	symbol	rpkm	matrix	partcount	geneexpress	cellcount	partid	pname
#tname	annotation	tfgenename	tfcount	rgenename	rcount	nprttype

print(dtfm)
print(dtfm['tname'])

#nan 삭제..
nprtdata = dtfm.dropna()
#column delete
#nprtdata =nprtdata.drop(['nprttype'],axis=1)
print("=============penguins==============")
print(nprtdata)
#tname	annotation	tissueid	tfgenename	rgenename	rcount
print(nprtdata.tname.value_counts())
print(nprtdata.annotation.value_counts())
print(nprtdata.tfgenename.value_counts())
print(nprtdata.rgenename.value_counts())


tnamedict ={}
npuniqtname= np.unique(nprtdata.tname)
print("======agenename[0]=========",npuniqtname[0])
for i in range(0,len(np.unique(nprtdata.tname))):
    tnamedict[npuniqtname[i]] = i

tfgenenamedict ={}
npuniqtfgenename= np.unique(nprtdata.tfgenename)
print("======agenename[0]=========",npuniqtfgenename[0])
for i in range(0,len(np.unique(nprtdata.tfgenename))):
    tfgenenamedict[npuniqtfgenename[i]] = i

rgenenamedict ={}
npuniqrgenename= np.unique(nprtdata.rgenename)
print("======agenename[0]=========",npuniqrgenename[0])
for i in range(0,len(np.unique(nprtdata.rgenename))):
    rgenenamedict[npuniqrgenename[i]] = i

annotationdict ={}
npuniqannotationname= np.unique(nprtdata.annotation)
print("======npuniqannotationname[0]=========",npuniqannotationname[0])
for i in range(0,len(np.unique(nprtdata.annotation))):
    annotationdict[npuniqannotationname[i]] = i


print(tnamedict)
print(tfgenenamedict,len(tfgenenamedict))
print(rgenenamedict,len(rgenenamedict))

tfgenenamedict1= tfgenenamedict.copy()
print("===========tfgenenamedict1=============")
#초기화..
tfgenenamedict1=tfgenenamedict1.fromkeys(tfgenenamedict1,0)
print(tfgenenamedict1.values())

rgenenamedict1=rgenenamedict.copy()
#dict.fromkeys(['X', 'Y', 'Z'], 0)
print("===========rgenenamedict1=============")
#초기화..
rgenenamedict1=rgenenamedict1.fromkeys(rgenenamedict1,0)
print(rgenenamedict1)

for key in rgenenamedict1:
    print(f"{key}: {rgenenamedict1[key]}")



print("=======nprtdata.columns=======")
print(nprtdata.shape)
print("=======nprtdata.columns=======")
print(nprtdata.columns)
print("=======nprtdata.describe()=======")
print(nprtdata.describe())
print("=======nprtdata.info()=======")
print(nprtdata.info())


tname_annotation = pd.DataFrame(
nprtdata.groupby(['tname', 'annotation','tfgenename','rgenename'])["rgenename"].count()).rename(columns={"rgenename": "Count"})

print("==========tname_annotation=======")
print(tname_annotation)


Body_df = nprtdata[
nprtdata['tname'] == 'Body'
].reset_index(drop=True)

print("======Body_df======")
print(Body_df)


#tname	annotation	tissueid	tfgenename	rgenename	rcount
print(Body_df.tname.value_counts())
print(Body_df.annotation.value_counts())
print(Body_df.tfgenename.value_counts())
print(Body_df.rgenename.value_counts())

#npgene sum(count)
#tname	annotation	tfgenename	tfcount	rgenename	rcount
rgenecountsum= Body_df.groupby('rgenename')['rcount'].sum()
motifnamesum= Body_df.groupby('tfgenename')['tfcount'].sum()


G = nx.Graph(day="Stackoverflow")

'''
for index, row in Body_df.iterrows():
    G.add_node(row['tfgenename'], group=row['annotation'], nodesize=row['rcount'])

for index, row in Body_df.iterrows():
    G.add_weighted_edges_from([(row['annotation'], row['tfgenename'], row['rcount'])])
'''

for index, row in Body_df.iterrows():
    #G.add_node(row['motifname'], group=row['npgenename'], nodesize=row['tfcount'])
    G.add_node(row['tfgenename'], group=row['rgenename'], nodesize=1000 )

for index, row in Body_df.iterrows():
    #G.add_node(row['motifname'], group=row['npgenename'], nodesize=row['tfcount'])
    G.add_node(row['rgenename'], nodesize=row['tfcount'])


for index, row in Body_df.iterrows():
    #G.add_weighted_edges_from([(row['tfgenename'], row['rgenename'], row['rcount'])])
    G.add_edges_from([(row['tfgenename'], row['rgenename'])])

    # tname	motifname	tfgenename	tfcount	np
    # genename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight

    #G.add_weighted_edges_from([(row['motifname'], row['npgenename'], row['npmotifregweight'])])
    #G.add_nodes_from([(row['motifname'], row['npgenename'], row['npmotifregweight'])])
    #G.add_edges_from([(row['motifname'], row['npgenename'])])




color_map = {1:'#f09494', 2:'#eebcbc', 3:'#72bbd0', 4:'#91f0a1', 5:'#629fff', 6:'#bcc2f2',
             7:'#eebcbc', 8:'#f1f0c0', 9:'#d2ffe7', 10:'#caf3a6', 11:'#ffdf55', 12:'#ef77aa',
             13:'#d6dcff', 14:'#d2f5f0'}

colors=[]
tfgenecount = Body_df.tfgenename.value_counts()
rgenecount= Body_df.rgenename.value_counts()

print("===Bodyrgenecountsum===")

print(rgenecountsum.index)

for i in rgenecountsum.index:
    print(i)
    print(rgenecountsum.get(i))

for i in motifnamesum.index:
    print(i)
    print(motifnamesum.get(i))

node_sizes=[]
for node in G:
    print("node,rgenecount.get(node)")
    print(node,rgenecount.get(node))
    if rgenecount.get(node) ==None:
        # colors.append([1, 1, 0.5])
        colors.append('#f09494')
        node_sizes.append(motifnamesum.get(node)/1000)
    else:
        # colors.append([1,0.5,0.5])
        colors.append('#caf3a6')
        node_sizes.append(rgenecountsum.get(node)/1000)


#colors = [color_map[G.node[node]['rgenename']] for node in G]
#sizes = [G.node[node]['nodesize']*10 for node in G]
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['font.family'] = 'Helvetica'  #'Helvetica'
plt.rcParams['font.size'] = 12

plt.figure(figsize=(13,13))
plt.title('Body', fontsize=20)
options = {
    'edge_color': '#FFDEA2',
    'width': 1,
    'with_labels': True,
    'font_weight': 'regular',
    'font_family': 'Helvetica',
    'font_size': 10,
}

f = plt.figure(1)
#ax = f.add_subplot(1,1,1)
ax = plt.gca()

ax.plot([0],[0],color='#FFDEA2',label='Edge')
ax.plot([0],[0],color='#f09494',label='TFs')
ax.plot([0],[0],color='#caf3a6',label='NP')


print("node size", node_sizes)
#legend_elements = [Line2D([0],[0],color='#FFDEA2',label='edge'),
#Line2D([0],[0],marker='o',markersize=10,color='#f09494',label='TF:tj_(+) size(7206)'),
#Line2D([0],[0],marker='o',markersize=10,color='#72bbd0',label='NPR:CrzR size(721)')

#nx.draw(G, node_color=colors, node_size=sizes,pos=nx.spring_layout(G, k=0.25, iterations=50), **options)
#nx.draw(G, node_color=colors,  pos=nx.spring_layout(G, k=1.0, iterations=6), **options)
#nx.draw_networkx(G, node_color=colors,  pos=nx.spring_layout(G, k=1, iterations=6), **options,ax=ax)
nx.draw_networkx(G, node_color=colors,  pos=nx.spring_layout(G, k=0.25, iterations=7), **options,ax=ax ,node_size=node_sizes)


f.set_facecolor('w')
plt.legend(fontsize=12)
f.tight_layout()
#plt.axis('off')
#nx.draw_networkx(G,pos=nx.spring_layout(G, k=0.25, iterations=10), node_color=colors,with_labels=True,ax=ax)
#ax = plt.gca()
#ax.collections[0].set_edgecolor("#555555")
#plt.legend(ncol=1,loc='right',labels=['transcription factor','edge line','neuro peptide'])

plt.show()




Head_df = nprtdata[
nprtdata['tname'] == 'Head'
].reset_index(drop=True)

print("======Head_df======")
print(Head_df)






G = nx.Graph(day="Stackoverflow")

'''
for index, row in Head_df.iterrows():
    G.add_node(row['rgenename'], group=row['annotation'], nodesize=row['rcount'])

for index, row in Head_df.iterrows():
    G.add_weighted_edges_from([(row['annotation'], row['rgenename'], row['rcount'])])
'''
#tname	annotation	tissueid	tfgenename	rgenename	rcount

for index, row in Head_df.iterrows():
    G.add_node(row['tfgenename'], group=row['rgenename'], nodesize=row['rcount'])

for index, row in Head_df.iterrows():
    G.add_weighted_edges_from([(row['tfgenename'], row['rgenename'], row['rcount'])])


color_map = {1:'#f09494', 2:'#eebcbc', 3:'#72bbd0', 4:'#91f0a1', 5:'#629fff', 6:'#bcc2f2',
             7:'#eebcbc', 8:'#f1f0c0', 9:'#d2ffe7', 10:'#caf3a6', 11:'#ffdf55', 12:'#ef77aa',
             13:'#d6dcff', 14:'#d2f5f0'}

for node in G:
    print("=====node======")
    print(node)


#if node_labels is not None:
#colors = [[1, 1, 0.5] if x else [1,0.5,0.5] for x in G ]
colors=[]
#print(colors)
tfgenecount = Head_df.tfgenename.value_counts()
rgenecount= Head_df.rgenename.value_counts()

#tname	annotation	tfgenename	tfcount	rgenename	rcount

rgenecountsum= Head_df.groupby('rgenename')['rcount'].sum()
motifnamesum= Head_df.groupby('tfgenename')['tfcount'].sum()

print("===rgenecountsum===")

print(rgenecountsum.index)

for i in rgenecountsum.index:
    print(i)
    print(rgenecountsum.get(i))

for i in motifnamesum.index:
    print(i)
    print(motifnamesum.get(i))


node_sizes = []

print("===Head_df.npmotifregoccur===")
print(Head_df.tfcount)
print(Head_df.rgenename)
#tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight
#tname	annotation	tfgenename	tfcount	rgenename	rcount


for node in G:
    print("======node=====")
    print(node,rgenecount.get(node))
    if rgenecount.get(node) ==None:
        #colors.append([1, 1, 0.5])
        colors.append('#f09494')
        node_sizes.append(motifnamesum.get(node)/1000)
    else:
        #colors.append([1,0.5,0.5])
        colors.append('#caf3a6')
        node_sizes.append(rgenecountsum.get(node)/1000)




#colors = [color_map[G.node[node]] for node in G]
#sizes = [G.node[node]['nodesize']*10 for node in G]



print(colors)
#colors = [color_map[G.node[node]] for node in G]
#sizes = [G.node[node]['nodesize']*10 for node in G]

plt.figure(figsize=(13,13))
plt.title('Head', fontsize=20)
options = {
    'edge_color': '#FFDEA2',
    'width': 1,
    'with_labels': True,
    'font_weight': 'regular',
    'font_family' : 'Helvetica',
    'font_size' : 10,
}

#nx.draw(G, node_color=colors, node_size=sizes,pos=nx.spring_layout(G, k=0.25, iterations=50), **options)
#nx.draw(G, pos=nx.spring_layout(G, k=0.15, iterations=5), **options)
#nx.draw(G, pos=nx.spring_layout(G, k=0.45, iterations=7), **options)

#if node_labels is not None:
#        node_color = [[1, 1, 0.5] if x else [1,0.5,0.5] for x in node_labels]

f = plt.figure(1)
#ax = f.add_subplot(1,1,1)
ax = plt.gca()
f.subplots_adjust(top=0.85)

plt.rcParams['lines.linewidth'] = 1
plt.rcParams['font.family'] = 'Helvetica'  #'Helvetica'
plt.rcParams['font.size'] = 12


ax.plot([0],[0],color='#FFDEA2',label='Edge')
ax.plot([0],[0],color='#f09494',label='TFs')
ax.plot([0],[0],color='#caf3a6',label='NP')


nx.draw_networkx(G, node_color=colors,  pos=nx.spring_layout(G, k=0.01 , iterations=5), **options , ax=ax , node_size=node_sizes)

f.set_facecolor('w')
plt.legend(fontsize=12)
f.tight_layout()


#ax = plt.gca()
#ax.collections[0].set_edgecolor("#555555")
#plt.legend(ncol=1,loc='right',labels=['transcription factor','edge line','neuro peptide'])
plt.show()

'''

G=nx.Graph()
G.add_node("kind1")
G.add_node("kind2")
G.add_node("kind3")
G.add_node("kind4")
G.add_node("kind5")
G.add_node("kind6")

# You were missing the position.
pos=nx.spring_layout(G)
val_map = {'kind1': 2,'kind2': 2,'kind3': 2,'kind4': 1,'kind5':4,'kind6': 3}
#I had this list for the name corresponding t the color but different from the node name
ColorLegend = {'Obsolete': 2,'Initialisation': 1,'Draft': 4,'Release': 3}
values = [val_map.get(node, 0) for node in G.nodes()]
# Color mapping
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=max(values))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# Using a figure to use it as a parameter when calling nx.draw_networkx
f = plt.figure(1)
ax = f.add_subplot(1,1,1)
for label in ColorLegend:
    ax.plot([0],[0],color=scalarMap.to_rgba(ColorLegend[label]),label=label)

# Just fixed the color map
nx.draw_networkx(G,pos, cmap = jet, vmin=0, vmax= max(values),node_color=values,with_labels=True,ax=ax)

# Setting it to how it was looking before.
plt.axis('off')
f.set_facecolor('w')

plt.legend()

f.tight_layout()
plt.show()

'''