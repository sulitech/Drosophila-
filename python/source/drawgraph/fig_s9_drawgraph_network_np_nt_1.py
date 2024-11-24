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
from matplotlib.lines import Line2D

import matplotlib.colors as colors
import matplotlib.cm as cmx

from matplotlib.pyplot import rc_context

#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns


#tf와 np 클러스터
#size up node display different

filepath="D:/project/python/bioinfo/data/loom"
#filepath="D:/project/python/bioinfo/data/h5ad"
filepath="D:/project/python/bioinfo/data/graph"
hg19="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/hg19/"
chgnew="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/chgnew"
pheromone="D:/project/python/bioinfo/data/graph/creatematrix/pheromone/"
creatematrixpath="D:/project/python/bioinfo/data/graph/creatematrix"
chemopath="D:/project/python/bioinfo/data/graph/chemo"
tfpath="D:/project/python/bioinfo/data/graph/nt"
savepath="D:/project/python/bioinfo/figures/nt/network/np"


filename="r_fca_biohub_testis_10x.loom"
filename="s_fca_biohub_antenna_10x.loom"
filename="s_fca_biohub_head_10x.loom"
filename="Chemosensory.txt"
filename="np_10_regulon_tissue_ann_tfgene_reggene_list.txt"
filename="np_nt_group.txt"   #np tissue 10

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

#For whole tissue

filedir= os.path.join(tfpath,filename)

print(filedir)

#pbarcodes = pd.DataFrame(os.path.join(chemopath,filename),index=0)
dtfm = pd.read_csv(os.path.join(tfpath,filename), delimiter='\t')
#tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight

print(dtfm)
print(dtfm['tname'])


#nan 삭제..
#nprtdata = dtfm.dropna()
nprtdata = dtfm
#column delete
#nprtdata =dtfm.drop(['nprttype'],axis=1)
print("=============penguins==============")
print(nprtdata)

print(nprtdata.tname.value_counts())
#print(nprtdata.motifname.value_counts())
#print(nprtdata.tfgenename.value_counts())
print(nprtdata.npgenename.value_counts())


tnamedict ={}
npuniqtname= np.unique(nprtdata.tname)
print("======agenename[0]=========",npuniqtname[0])
for i in range(0,len(np.unique(nprtdata.tname))):
    tnamedict[npuniqtname[i]] = i

tfgenenamedict ={}


#tf->nt
npuniqtfgenename= np.unique(nprtdata.nt)
print("======agenename[0]=========",npuniqtfgenename[0])
for i in range(0,len(np.unique(nprtdata.nt))):
    tfgenenamedict[npuniqtfgenename[i]] = i

rgenenamedict ={}
npuniqrgenename= np.unique(nprtdata.npgenename)
print("======agenename[0]=========",npuniqrgenename[0])
for i in range(0,len(np.unique(nprtdata.npgenename))):
    rgenenamedict[npuniqrgenename[i]] = i


'''
annotationdict ={}
npuniqannotationname= np.unique(nprtdata.annotation)
print("======npuniqannotationname[0]=========",npuniqannotationname[0])
for i in range(0,len(np.unique(nprtdata.annotation))):
    annotationdict[npuniqannotationname[i]] = i
'''
#tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight


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

'''
tname_annotation = pd.DataFrame(
nprtdata.groupby(['tname', 'annotation','tfgenename','rgenename'])["rgenename"].count()).rename(columns={"rgenename": "Count"})
print("==========tname_annotation=======")
print(tname_annotation)
'''

ltname=['Testis','Antenna','Body wall','Fat body','Gut','Haltere','Heart','Leg','Male reproductive glands','Malpighian tubule',
        'Oenocyte','Ovary','Proboscis','Trachea','Wing','Head','Body']


def drawGraph(ltname):


    for j in range(17):

        print("===i===",j)
        stname = ltname[j]

        Body_df = nprtdata[
            nprtdata['tname'] == stname
            ].reset_index(drop=True)

        print("======Body_df======1 ")
        print(Body_df)


        # tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight

        print(Body_df.tname.value_counts())
        #print(Body_df.motifname.value_counts())
        #print(Body_df.tfgenename.value_counts())
        print(Body_df.npgenename.value_counts())

        # npgene sum(count)
        rgenecountsum = Body_df.groupby('npgenename')['npexpress'].sum()
        motifnamesum = Body_df.groupby('nt')['ntexpress'].sum()


        G = nx.Graph(day="Stackoverflow")

        '''
        for index, row in Body_df.iterrows():
            G.add_node(row['tfgenename'], group=row['annotation'], nodesize=row['rcount'])

        for index, row in Body_df.iterrows():
            G.add_weighted_edges_from([(row['annotation'], row['tfgenename'], row['rcount'])])
        '''
        # tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight
        for index, row in Body_df.iterrows():
            # G.add_weighted_edges_from([(row['motifname'], row['npgenename'], row['npmotifregweight']*100)])
            G.add_node(row['npgenename'], nodesize=row['npexpress'])


        for index, row in Body_df.iterrows():
            # G.add_node(row['motifname'], group=row['npgenename'], nodesize=row['tfcount'])
            G.add_node(row['nt'], nodesize=10)


        for index, row in Body_df.iterrows():
            G.add_edges_from([(row['nt'], row['npgenename'])])
            # G.add_nodes_from([(row['motifname'], row['npgenename'], row['npmotifregweight'])])
            # G.add_edges_from([(row['motifname'], row['npgenename'])])

        color_map = {1: '#f09494', 2: '#eebcbc', 3: '#72bbd0', 4: '#91f0a1', 5: '#629fff', 6: '#bcc2f2',
                     7: '#eebcbc', 8: '#f1f0c0', 9: '#d2ffe7', 10: '#caf3a6', 11: '#ffdf55', 12: '#ef77aa',
                     13: '#d6dcff', 14: '#d2f5f0'}
        # tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight

        colors = []
        tfgenecount = Body_df.nt.value_counts()
        rgenecount = Body_df.npgenename.value_counts()

        print("===Bodyrgenecountsum===")

        print(rgenecountsum.index)

        for i in rgenecountsum.index:
            print(i)
            print(rgenecountsum.get(i))

        for i in motifnamesum.index:
            print(i)
            print(motifnamesum.get(i))

        node_sizes = []

        for node in G:
            print(node, rgenecount.get(node))
            if rgenecount.get(node) == None:
                # colors.append([1, 1, 0.5])

                if node == motifnamesum.index[0]:
                    colors.append('#7e4e7e')
                else:
                    colors.append('#d8bfd8')

                node_sizes.append(motifnamesum.get(node) / 10)

            else:
                # colors.append([1,0.5,0.5])

                if node == rgenecountsum.index[0]:
                    colors.append('#3d700f')
                else:
                    colors.append('#caf3a6')

                node_sizes.append(rgenecountsum.get(node)/10)

        # colors = [color_map[G.node[node]['rgenename']] for node in G]
        # sizes = [G.node[node]['nodesize']*10 for node in G]
        plt.rcParams['lines.linewidth'] = 1
        plt.rcParams['font.family'] = 'Helvetica'  # 'Helvetica'
        plt.rcParams['font.size'] = 12

        plt.figure(figsize=(13, 13))
        plt.title( stname , fontsize=20)

        options = {
            'edge_color': '#FFDEA2',
            'width': 1,
            'with_labels': True,
            'font_weight': 'regular',
            'font_family': 'Helvetica',
            'font_size': 10,
        }

        f = plt.figure(1)
        # ax = f.add_subplot(1,1,1)
        ax = plt.gca()

        #print('====rgenecountsum.index(0)====')
        #print(rgenecountsum.index[0])
        #print(rgenecountsum.get(0))

        #print(motifnamesum.index[0])
        #print(motifnamesum.get(0))
        #('He is %s %s and is %i years old.') % (first_name, last_name, age)


        legend_elements = [Line2D([0], [0], color='#FFDEA2', label='Edge'),
                           Line2D([0], [0], marker='o', markersize=10, color='#d8bfd8', label='NT:'+motifnamesum.index[0]+' size('+str(motifnamesum.get(0))+')'),
                           Line2D([0], [0], marker='o', markersize=10, color='#caf3a6', label='NP:'+rgenecountsum.index[0]+' size('+str(rgenecountsum.get(0))+')')
                           ]
        # nx.draw(G, node_color=colors, node_size=sizes,pos=nx.spring_layout(G, k=0.25, iterations=50), **options)
        # nx.draw(G, node_color=colors,  pos=nx.spring_layout(G, k=1.0, iterations=6), **options)
        nx.draw_networkx(G, node_color=colors, pos=nx.spring_layout(G, k=0.5, iterations=50), **options, ax=ax,
                         node_size=node_sizes)

        f.set_facecolor('w')
        plt.legend(handles=legend_elements, fontsize=12, loc='upper left')
        f.tight_layout()

        savepath = "D:/project/python/bioinfo/figures/nt/network/np"
        #savefile = "figure5_nt_np_network_"+ ltname[j]+".png"
        savefile = "figure5_nt_np_network_"+ ltname[j]+".pdf"
        plt.savefig(os.path.join(savepath, savefile))

        plt.show()

        print(ltname[j])


#each tissue shows tf, np
drawGraph(ltname)



