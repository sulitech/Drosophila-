import loompy
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import scipy

#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns

# this page contains handling loom file

filepath="D:/project/python/bioinfo/data"
filename="r_fca_biohub_testis_10x.loom"

partid = 8
#filename="s_fca_biohub_head_10x.h5ad"

if partid==3:
    clockpath="D:/project/python/bioinfo/data/graph/clock/body/"
elif partid==8:
    clockpath = "D:/project/python/bioinfo/data/graph/clock/head/"



filedir= os.path.join(filepath,filename)

print(filedir)

#with loompy.connect(filedir) as ds:
#    print(ds)

#adata = anndata.AnnData

#print("===========adata=================")
#print(adata)
#adata.X =

ds = loompy.connect(filedir)
print("===========ds=================")
#12158 6527
print(ds.shape)
print(ds)


nparr = np.array(ds[:, :], dtype='int')
print("=============ds nparr======================")
print(nparr)
anndataX= scipy.sparse.csc_matrix(nparr)

#scipy.io.mmwrite(os.path.join(pheromone,mfilename), narry)

# load sparse matrix:
#X = io.mmread("counts.mtx")
# create anndata object
#adata = anndata.AnnData(
#    X=X.transpose().tocsr()
#)


print("=============h5ad anndataadata.shape======================")

print(anndataX.shape)

#0,0 are data, if change to mtx  1,1 

totals = ds.map([np.sum], axis=1)[0]  # Calculate the total molecule count for each cell
print("total==>",totals)
cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
print("cells==>",len(cells))


print("=====col_attrs==============")
print(ds.col_attrs)
print("=====row_attrs==============")
print(ds.row_attrs)
#totals = ds.map([np.sum], axis=1)[0]
#cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
print("========tissue===========")
print(ds.col_attrs['tissue'])

#print(ds.row_attrs)
print("=======gene==============")
gene= np.array(ds.row_attrs['Gene'])

print(gene)
print(len(gene))
print(gene[0].decode('utf8'))

print("===========barcode==============")
print(ds.col_attrs['Barcode'])
print(len(ds.col_attrs['Barcode']))

print("===========annotation__ontology_id==============")
print(ds.col_attrs['annotation__ontology_id'])
print(len(ds.col_attrs['annotation__ontology_id']))

print("===========cell id==============")
print(ds.col_attrs['CellID'])
print(len(ds.col_attrs['CellID']))
print("===========annontation===========")
print(ds.col_attrs['annotation'])
print(len(ds.col_attrs['annotation']))  

print("===========MotifRegulons================")
print(ds.row_attrs['MotifRegulons'])
print(len(ds.row_attrs['MotifRegulons'])) .

print(ds.row_attrs['MotifRegulons'])
print(len(ds.row_attrs['MotifRegulons']))  

print("===========MotifRegulons= Abd-B_(+)-motif MotifRegulons===============")
print(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif'])  
print(ds.row_attrs['MotifRegulons'].dtype)  #314

print("===========MotifRegulons dtype==============")
#print(ds.row_attrs['MotifRegulons'].dtype)
#tf = np.array(ds.row_attrs['MotifRegulons'].dtype)

tf = [x for x,y in ds.row_attrs['MotifRegulons'].dtype.fields.items()]
print(len(tf))
print(tf[0])  #motif full name
print(tf[0].index('_('))
print(tf[0][:tf[0].index('_(')])  #motif gene name
#'MotifRegulonGeneWeights
print(ds.row_attrs['MotifRegulons'][tf[0]]) #matrix

#순서
print(ds.row_attrs['MotifRegulons'][tf[0]][20]) #matrix
print(ds.row_attrs['MotifRegulonGeneWeights'][tf[0]][20]) 
print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'][12094])

motifregeulons = np.transpose(np.nonzero(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif']))
print(motifregeulons)



motifregeulons1=np.transpose(np.nonzero(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif']))
motifregeulonsw= np.transpose(np.nonzero(ds.row_attrs['MotifRegulonGeneWeights']['Abd-B_(+)-motif']))
#print(motifregeulons1)
print(len(motifregeulons1))
print(len(motifregeulonsw))
#print(motifregeulons[20])
MotifRegulonGeneOccurrences = np.transpose(np.nonzero(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif']))
print(len(MotifRegulonGeneOccurrences))
print(ds.row_attrs['MotifRegulonGeneOccurrences'].shape)  #9076
print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'])

print(motifregeulons1[0][0])
print(motifregeulonsw[0][0])
print(MotifRegulonGeneOccurrences[0][0])


print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'][12094])
print(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif'][12094]) 
print(ds.row_attrs['MotifRegulonGeneWeights']['Abd-B_(+)-motif'][12094]) 

print(ds.col_attrs['MotifRegulonsAUC'].dtype)
print(ds.col_attrs['MotifRegulonsAUC'].shape)

tfauc = [x for x,y in ds.col_attrs['MotifRegulonsAUC'].dtype.fields.items()]
print(len(tfauc)) #314
MotifRegulontfauc = np.transpose(np.nonzero(ds.col_attrs['MotifRegulonsAUC']['Abd-B_(+)-motif']))
print(ds.col_attrs['MotifRegulonsAUC']['Abd-B_(+)-motif'][0])
#6517
print(len(MotifRegulontfauc))
print(len(ds.col_attrs['MotifRegulonsAUC']['Abd-B_(+)-motif']))
print(len(ds.col_attrs['MotifRegulonsAUC']))
#6527


#adata =anndata.AnnData(anndataX.T,obs=ds.col_attrs , var = ds.row_attr)

#obs_names = "cell_names",
#var_names = "gene_names",


print("================adata var===============")
#print(adata.var)

#Abd-B_(+)-motif
'''
get_regulons(
  loom,
  column.attr.name = "MotifRegulons",
  tf.as.name = TRUE,
  tf.sep = "_"
)
# Regulons:
regulons_incidMat <- get_regulons(loom) # as incid matrix
regulons_motif <- regulonsToGeneLists(regulons_incidMat) # convert to list
regulons_ChIP <- regulonsToGeneLists(get_regulons(loom, attrName = "TrackRegulons"))
# Regulon AUC and thresholds
regulonsAUC <- get_regulonsAuc(loom)
regulonsAucThresholds <- get_regulonThresholds(loom)
# Embeddings (tsne/umap)
embeddings <- get_embeddings(loom)
'''
print("================adata obs===============")
#print(adata.obs)

ds.close()
