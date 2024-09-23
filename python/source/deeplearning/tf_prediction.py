from gensim.models import Word2Vec
import numpy as np
from Cython.Build import cythonize
import logging
import os
import sys

print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))  #D:\data\slb\project\python
basedir = os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))) #D:\data\slb\project\python
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
sys.path.append(os.path.dirname(os.path.dirname(__file__))) #load_info

import sysinfo
import platform
import pymysql



dirpath = "D:/project/python/bioinfo/data/graph/creatematrix/nprt"
modelpath = "D:/project/python/bioinfo/model"
matrixdir = "D:/project/python/bioinfo/data/graph/creatematrix/nprt/matrix/"



logging.basicConfig(level=logging.INFO)

patterns = []




def load_sequence(from_path):
    with open(from_path) as fp:
        [patterns.append(line.strip().split("\t")) for line in fp]
        #[patterns.append(line.strip().split(",")) for line in fp]


def main():
    #train_path = "D:/data/slb/project/python/product_recom/input/tr_set_product_recom.txt"

    trainfilename = "traingene"  #annotatin없는것
    trainfilename = "npr_tf_tissue_training_tab"
    #trainfilename = "training_gene_np_npr_nt"
    #trainfilename ="tissue_ann_traingene"  #annotation있는것
    #trainfile = tissue_ann_traingene.txt

    train_path = os.path.join(modelpath,trainfilename+".txt")
    model_data_path=os.path.join(modelpath, trainfilename+".model")

    # curs = conn.cursor(pymysql.cursors.DictCursor)
    # SQL문 실행
    #sql = " truncate table product_recom"
    #curs.execute(sql)

    # split patterns to train_patterns and test_patterns


    print("model_data_path====>",model_data_path)

    if os.path.isfile(model_data_path):
        model = Word2Vec.load(model_data_path)

    #model = Word2Vec.load("frequent_model.model")

    # Test
    #test_size = float(len(test_patterns))
    hit = 0.0
    findcount =0

    #last_item = ['MsR1']
    last_item=['T_Atf-2','T_Atf3','T_Camta','T_Cnot4','T_Crtc','T_E2f1','T_Mrtf','T_Mst84B','T_T48','T_TfIIB',	'T_Tfb4','T_pre-rRNA:CR45845']
    last_item=['T_Atf3','T_Cnot4','T_E2f1','T_MTF-1','T_Mst36Fa']
    last_item = ['Ms', 'MsR2']
    #last_item=['Trissin','TrissinR']
    last_item=['Acp26Aa','SPR']
    last_item=['Akh','AkhR']

    last_item=['AstA','AstA-R1','AstA-R2']
    last_item=['slp1',	'Tet'	,'tj'	,'zfh2']
    last_item=['Atf3','CrebA','Hnf4','hth','Rel','sima','SREBP','srp','svp']

    last_item = ['svp' ,'onecut' , 'dimm' ]
    last_item = ['slou', 'Ptx1', 'Lmx1a']
    last_item = ['hth', 'pan', 'Rx', 'E2f1', 'pho']




    #last_item = ['Atf3','CrebA','Hnf4','hth','Rel','sima','SREBP']
    #last_item=['bru3',	'Cf2',	'CG16779',	'CG44247',	'CHES-1-like',	'crp',	'ey']
    #last_item=['acj6',	'Atac3',	'bigmax'	,'bru3',	'Cf2'	,'CG16779'	,'CG44247',	'CG9727']
    #last_item=['foxo']
    #last_item = ['Ilp6', 'InR']
    #last_item = ['CrebA',	'Hnf4',	'hth','Rel',	'sima',	'SREBP',	'srp',	'svp']
    #last_item = ['Camta','Cnot4','Mrtf','ftz-f1']


    try:
        #prediction = model.wv.most_similar(positive=last_item, topn=20)
        #prediction = model.wv.most_similar(positive=last_item, topn=30)
        prediction = model.wv.most_similar(negative=last_item, topn=50)

        # Check if the item that we have removed from the test, last_item, is among
        # the predicted ones.
        #print("prediction===>", prediction)
        for predicted_item, score in prediction:
            print("predicted_item", predicted_item)
            print("score", score)
            if predicted_item.find('T_')>=0 and score >=0.5:
                #print("last_item==>", last_item)
                print("predicted_item", predicted_item)
                print("score", score)

    except Exception as ex:
        print("erro=>",ex)


if __name__ == '__main__':
    main()
