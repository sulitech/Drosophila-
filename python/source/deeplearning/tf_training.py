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
dirpath="D:/project/python/bioinfo/data/graph/tf"
modelpath = "D:/project/python/bioinfo/model"
matrixdir = "D:/project/python/bioinfo/data/graph/creatematrix/nprt/matrix/"


logging.basicConfig(level=logging.INFO)

patterns = []

def load_sequence(from_path):
    icnt = 0
    with open(from_path) as fp:
        for line in fp:
            lsplits = line.strip().split("\t")
            if icnt !=0:
                lists=[]
                for lsplit in lsplits:
                    print("lsplit==>", lsplits[0],lsplit)
                    lists.append(lsplit)
                    lists.append(lsplits[0])
                patterns.append(lists)
            icnt=icnt+1



'''
        [patterns.append(line.strip().split("\t")) for line in fp]
        #[patterns.append(line.strip().split(",")) for line in fp]
    text = f.read()
    words = text[:10000].split()
    stop_words = set(stopwords.words("english"))
    words_filter = [w for w in words if not w in stop_words]
'''


def main():
    trainfilename = "traingene"  #annotatin없는것
    trainfilename = "training_gene_np_npr_tf"
    trainfilename = "training_gene_np_npr_nt"
    trainfilename="npr_tf_tissue_training_tab"
    #trainfilename ="tissue_ann_traingene"  #annotation있는것
    #trainfile = tissue_ann_traingene.txt

    train_path = os.path.join(dirpath,trainfilename+".txt")
    load_sequence(train_path)
    # split patterns to train_patterns and test_patterns

    print(patterns)


    #train_patterns = np.random.choice(patterns, np.floor(len(patterns) * 0.8))
    #test_patterns = np.random.choice(patterns, np.floor(len(patterns) * 0.2))
    train_size = int(len(patterns) * 0.7)
    #train_size = int(len(patterns))
    test_size = int(len(patterns) - train_size)

    train_size = int(len(patterns))
    #train_patterns = np.random.choice(patterns, train_size)
    test_patterns = np.random.choice(patterns, test_size)

    #train_patterns = int(len(patterns) * 0.7)
    #test_patterns = len(patterns) - train_size
    #train, test = datas[:, :d_length], datas[train_size:len(datas), :d_length]

    # train,test=datas[0:train_size,:29],datas[train_size:len(datas),:29]

    # Word vector representation learning
    #model = Word2Vec(patterns, window=20, min_count=2 , vector_size=35, workers=4, sample=1e-4 , hs=1, sg=0, epochs=1000)

    model = Word2Vec(min_count=2,
                         window=100,
                         sample=6e-5,
                         alpha=0.03,
                         min_alpha=0.0007,
                         negative=50,
                         vector_size=25,
                         workers=4,
                     sg=0)

    print("=======model.build_vocab(patterns, progress_per=10000)============")
    model.build_vocab(patterns, progress_per=10000)
    # window=100, vector_size=35,
    #min_count 최소단어 , window->앞뒤로 고려하는 단어수 140으로 했음..  vector_size 단어수*0.25 ,

    model.train(patterns ,total_examples=model.corpus_count, epochs=1000, report_delay=1)

    model.save(os.path.join(modelpath, trainfilename+".model"))

    # Test
    test_size = float(len(test_patterns))
    hit = 0.0

    print(test_patterns)

    #prediction = model.wv.most_similar(positive=['Gpb5'], topn=10)

    #print(prediction)

    for current_pattern in test_patterns:
        if len(current_pattern) < 2:
            test_size -= 1.0
            continue
        # Reduce the current pattern in the test set by removing the last item
        last_item = current_pattern.pop()

        print(last_item)
        # Keep those items in the reduced current pattern, which are also in the models vocabulary
        items = [it for it in current_pattern if it in model.wv.vocab]

        if len(items) <= 2:
            test_size -= 1.0
            continue

        # Predict the most similar items to items
        print("items===>",items)

        prediction = model.most_similar(positive=items , topn=3)

        # Check if the item that we have removed from the test, last_item, is among
        # the predicted ones.
        print("prediction===>", prediction)
        for predicted_item, score in prediction:
            if predicted_item == last_item:
                hit += 1.0
                #print("last_item==>", last_item)
                #print("prediction==>", prediction)
            else:
                print("last_item==>", last_item)
                print("predicted_item", predicted_item)
                print("score", score)



    print('Accuracy like measure: {}'.format(hit / test_size))


if __name__ == '__main__':
    main()
