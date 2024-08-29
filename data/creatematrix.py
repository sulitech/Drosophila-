import loompy
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import pymysql
import sysinfo

from matplotlib.pyplot import rc_context

#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns

filepath="D:/project/python/bioinfo/data/loom"
#filepath="D:/project/python/bioinfo/data/h5ad"
filepath="D:/project/python/bioinfo/data/graph/creatematrix"
hg19="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/hg19/"
pheromone="D:/project/python/bioinfo/data/graph/creatematrix/pheromone/"
chgnewdir="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/chgnew/"
filename="r_fca_biohub_testis_10x.loom"
filename="s_fca_biohub_antenna_10x.loom"
filename="s_fca_biohub_head_10x.loom"
mfilename="matrix.mtx"
gfilename="genes.tsv"
bfilename="barcodes.tsv"
#filename="s_fca_biohub_head_10x.h5ad"



conn = pymysql.connect(host=sysinfo.dburl, user=sysinfo.dbuser, password=sysinfo.dbpw,
                       db=sysinfo.database, charset='utf8')
curs = conn.cursor()


sql = " select concat(left(b.cellname,14),'-1') from mcell a , cell b "
sql = sql+" where a.cellid = b.cellid order by a.mcellid"

curs.execute(sql)
rows = curs.fetchall()

f = open(os.path.join(pheromone,bfilename),"w+")
for data in rows:
    f.write(data[0])
    f.write("\r")



sql = " select b.genename ,  b.symbol from mgene a ,  gene b "
sql = sql+" where a.geneid = b.geneid "
sql = sql+" order by a.geneid "



curs.execute(sql)
rows = curs.fetchall()

f = open(os.path.join(pheromone,gfilename),"w+")
for data in rows:
    f.write("{}\t{}\r".format(data[1], data[0]))


'''
sql = " select a.geneid, a.cellid,  a.geneexpress from part_cell_gene a "
sql = sql+" inner join gene b "
sql = sql+" on a.geneid = b.geneid "
sql = sql+" where a.tissueid =1 and a.cellid <=100 "
sql = sql+" order by a.geneid desc, a.cellid asc "
'''

sql = " select c.mgeneid , b.mcellid  ,a.geneexpression from  expgene a , "
sql = sql+" mcell b , mgene c "
sql = sql+" where a.cellid = b.cellid "
sql = sql+" and a.geneid = c.geneid "


#sql = sql+" order by rand() "
#sql = sql+" limit 1, 280000 "

curs.execute(sql)
rows = curs.fetchall()

f = open(os.path.join(pheromone,mfilename),"w+")

narry = np.array(rows)
print(narry.shape)
mxgene = np.max(narry[:,0])
mxcell = np.max(narry[:,1])
lencount=len(narry)

print(mxgene,mxcell,lencount)

f.write("%%MatrixMarket matrix coordinate real general\n")
f.write("%\n")
f.write(str(int(mxgene)) + " " + str(int(mxcell)) + " " + str(lencount) + "\n")
for data in rows:
    f.write("{} {} {}\n".format(data[0], data[1], data[2]))













f.close()
curs.close()
conn.close()
