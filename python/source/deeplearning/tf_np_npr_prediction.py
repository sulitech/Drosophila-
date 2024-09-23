import logging
import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


import keras
from keras.datasets import cifar10
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
# from keras.optimizers import SGD
from keras.constraints import maxnorm
from keras.models import load_model
from sklearn.metrics import roc_curve

import tensorflow as tf

batch_size = 28
num_classes = 1001
epochs = 1000
data_augmentation = True
num_predictions = 1001
save_dir = os.path.join(os.getcwd(), 'saved_models')
imageSize = 28

print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))  # D:\data\slb\project\python
basedir = os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))  # D:\data\slb\project\python
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
sys.path.append(os.path.dirname(os.path.dirname(__file__)))  # load_info

import sysinfo
import platform
import pymysql

tfpath = "D:/project/python/bioinfo/data/graph/tf"
modelpath = "D:/project/python/bioinfo/model"

filename = "r_fca_biohub_testis_10x.loom"
filename = "s_fca_biohub_antenna_10x.loom"
filename = "s_fca_biohub_head_10x.loom"
filename = "Chemosensory.txt"

npfilename = "np_tf_tissue_training.txt"
npnprfilename = "np_npr_tf_tissue_training.txt"
nprfilename = "npr_tf_tissue_training.txt"
savetraining = "save_np_npr_tf_tissue_training.pkl"

modeltname = "tname_training_tf" + str(epochs) + '.h5'
modelrgenename = "rgenename_training_tf" + str(epochs) + '.h5'

logging.basicConfig(level=logging.INFO)

patterns = []


def load_sequence(from_path):
    icnt = 0
    with open(from_path) as fp:
        for line in fp:
            lsplits = line.strip().split("\t")
            if icnt != 0:
                lists = []
                for lsplit in lsplits:
                    print("lsplit==>", lsplits[0], lsplit)
                    lists.append(lsplit)
                    lists.append(lsplits[0])
                patterns.append(lists)
            icnt = icnt + 1


filedir = os.path.join(tfpath, npfilename)

dtnp = pd.read_csv(os.path.join(tfpath, npfilename), delimiter='\t')
print(len(dtnp))

dtnpr = pd.read_csv(os.path.join(tfpath, nprfilename), delimiter='\t')
print(len(dtnpr))

dtnpnpr = pd.read_csv(os.path.join(tfpath, npnprfilename), delimiter='\t')
print(len(dtnpnpr))

# dtfm['np_npr'] = dtfm['npgenename'].map(str) + '_' + dtfm['nprgenename'].map(str)
# nan 삭제..
# nprtdata = dtfm.dropna()
# column delete
# nprtdata =nprtdata.drop(['nprttype'],axis=1)


dtnpnpr_np = dtnpnpr.copy()
dtnpnpr_np = dtnpnpr_np.drop(['rgenename'], axis=1)
dtnpnpr_np.rename(columns={'pgenename': 'rgenename'}, inplace=True)

print(dtnpnpr_np)

dtnpnpr_npr = dtnpnpr.copy()
dtnpnpr_npr = dtnpnpr_npr.drop(['pgenename'], axis=1)

print(dtnpnpr_npr)

dtnprtrain = pd.concat([dtnp, dtnpr, dtnpnpr_np, dtnpnpr_npr])
dtnprtrain.reset_index(inplace=True)

tnamedict = {}
npuniqtname = np.unique(dtnprtrain.tname)
print("======agenename[0]=========", npuniqtname[0])
for i in range(0, len(np.unique(dtnprtrain.tname))):
    tnamedict[npuniqtname[i]] = i

tfgenenamedict = {}
npuniqtfgenename = np.unique(dtnprtrain.tgenename)
print("======agenename[0]=========", npuniqtfgenename[0])
for i in range(0, len(np.unique(dtnprtrain.tgenename))):
    tfgenenamedict[npuniqtfgenename[i]] = i

rgenenamedict = {}
npuniqrgenename = np.unique(dtnprtrain.rgenename)
print("======agenename[0]=========", npuniqrgenename[0])
for i in range(0, len(np.unique(dtnprtrain.rgenename))):
    rgenenamedict[npuniqrgenename[i]] = i

print(tnamedict)
print(tfgenenamedict)
print(rgenenamedict)

print(dtnprtrain.tgenename)

tfgenenamedict = {}
alist = []
for text in dtnprtrain.tgenename:
    strings = text.split(',')
    for i in strings:
        alist.append(i)
# print(alist)


npuniqtfgenename = np.unique(alist)
print(npuniqtfgenename)
for i in range(0, len(npuniqtfgenename)):
    tfgenenamedict[npuniqtfgenename[i]] = i
print(tfgenenamedict)

print(len(dtnprtrain.tname))
print(len(dtnprtrain.tgenename))
print(len(dtnprtrain.rgenename))

print(len(dtnprtrain.columns))
print(len(dtnprtrain.values))

print("==========dtnprtrain.rgenename[0])========")
print(dtnprtrain)
print(dtnprtrain.loc[1]['rgenename'])
print(dtnprtrain.loc[1]['tname'])
print(dtnprtrain.loc[1]['tgenename'])
print(tnamedict[dtnprtrain.loc[1]['tname']])
print(rgenenamedict[dtnprtrain.loc[1]['rgenename']])

# print(dtnprtrain.rgenename[1])
# print(dtnprtrain.tname[0])


training = []
icount = 0

'''
for text in dtnprtrain.tgenename:
    strings = text.split(',')
    tfgenenamedict1 = tfgenenamedict.copy()
    # dict.fromkeys(['X', 'Y', 'Z'], 0)
    print("===========rgenenamedict1=============")
    # 초기화..
    tfgenenamedict1 = tfgenenamedict1.fromkeys(tfgenenamedict1, 0)
    print(tfgenenamedict1)
    print(strings)
    for i in strings:
        print(i)
        tfgenenamedict1[i]=1
        alist.append(i)

    print(tfgenenamedict1)
    list1 = list(tfgenenamedict1.values())

    list1.append(rgenenamedict[dtnprtrain.loc[icount]['rgenename']]) #append rgenename
    list1.append(tnamedict[dtnprtrain.loc[1]['tname']]) #append tname
    #training.append([list(tfgenenamedict1.values()),1,1])
    training.append(list1)
    icount=icount+1

    #print(training)
print("===========training=============")
print(training)
#arrtraining= np.array(training)
#print(arrtraining)
with open(os.path.join(tfpath,savetraining), 'wb') as save_training_file:
    # Step 3
    pickle.dump(training, save_training_file)

'''

with open(os.path.join(tfpath, savetraining), 'rb') as save_training_file:
    # Step 3
    training = pickle.load(save_training_file)

    # After config_dictionary is read from file
    print(training)

print(tnamedict)
print(tfgenenamedict)
print(rgenenamedict)

#tnamegetdict = dict((y,x) for x,y in tnamedict.iteritems())
#rgenenamegetdict = dict((y,x) for x,y in rgenenamedict.iteritems())

tnamegetdict = {y:x for x,y in tnamedict.items()}
rgenenamegetdict = {y:x for x,y in rgenenamedict.items()}

print("=========tnamegetdict====================")
print(tnamegetdict)
print("=========rgenenamegetdict====================")
print(rgenenamegetdict)



arrnp = np.array(training)

print(arrnp)
print(arrnp.shape)
print(len(arrnp))  # row
print(len(arrnp[0]))  # column

arrtraining = arrnp[:, 0:450]
print(len(arrtraining))
print(len(arrtraining[0]))

labelrgenename = arrnp[:, 450]
labeltname = arrnp[:, 451]

print(labelrgenename.shape)
print(labeltname.shape)

# x_train, y_train , labelMap =load_plateimage.getDataLabel(input)
# x_test, y_test,labelMap =load_plateimage.getDataLabel(inputt)


# Convert class vectors to binary class matrices.
# y_train = keras.utils.to_categorical(y_train, num_classes)
# y_test = keras.utils.to_categorical(y_test, num_classes)


print("===============num_classes============")
print(num_classes)

print(arrtraining.shape[1:])

label = labelrgenename
label = labeltname

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(arrtraining, labeltname, test_size=0.2,
                                                    random_state=1)
print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

#print(X_train)  # train value
#print(X_test)
#print(y_train)  # train label
#print(y_test)


X_train, X_test, y_train, rgene_test = train_test_split(arrtraining, labelrgenename, test_size=0.2,
                                                    random_state=1)
#print(X_train.shape, X_test.shape, y_train.shape, rgene_test.shape)

#print(X_train)  # train value
#print(X_test)
#print(y_train)  # train label
#print(rgene_test)

num_classes = len(rgenenamedict)
num_classes = len(tnamedict)

'''
model.fit( arrtraining , labelrgenename,
          batch_size=batch_size,
          epochs=epochs,
          shuffle=True)
'''
# label = labeltname
label = labelrgenename
label = labeltname


# Save model and weights
modeltname = "tname_training_tf" + str(epochs) + '.h5'
modelrgenename = "rgenename_training_tf" + str(epochs) + '.h5'


tnamemodel = load_model(os.path.join(modelpath, modeltname))

ty = tnamemodel.predict(X_test, verbose=0)
loss, accuracy = tnamemodel.evaluate(X_test, y_test)
print("tf Accuracy = {:.2f}".format(accuracy))
predictedtname = ty.argmax(axis=-1)
print(predictedtname)

rgenenamemodel = load_model(os.path.join(modelpath, modelrgenename))
ry=rgenenamemodel.predict(X_test, verbose=0)

print("=========X_test=============")
print(X_test)
predictedrgene = ry.argmax(axis=-1)
print(predictedrgene)
loss, accuracy = rgenenamemodel.evaluate(X_test, rgene_test)
print("rgene=====Accuracy=======")
print("Accuracy = {:.2f}".format(accuracy))


'''
precision    recall 
history = rgenenamemodel.fit(X_train, y_train,validation_split = 0.1, epochs=50, batch_size=4)
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'val'], loc='upper left')
plt.show()
'''

import numpy as np
from sklearn import datasets, metrics, model_selection, svm
from sklearn.metrics import classification_report,confusion_matrix ,plot_confusion_matrix



#y_test_class = np.argmax(y_test,axis=1)
#print(y_test_class)
#y_pred_class = np.argmax(ty,axis=1)
#print(y_pred_class)
print("======classification_report(y_test,predictedtname)=======")
print(classification_report(y_test,predictedtname))
print(confusion_matrix(y_test,predictedtname))
print("======classification_report(rgene_test,predictedrgene)=======")
print(classification_report(rgene_test,predictedrgene))
print("======confusion_matrix(rgene_test,predictedrgene)=======")
print(confusion_matrix(rgene_test,predictedrgene))



#X_train, X_test, y_train, rgene_test



txt="Atf3,bin,BtbVII,Cf2,CHES-1-like,cnc,CrebA,CrebB,Dif,foxo,ham,Hand,hb,Hnf4,jim,kay,Mef2,Myc,nej,REPTOR-BP,sim,so,SREBP,srp,svp,Trf2"
txt="Atf3,bin,BtbVII,Cf2,CHES-1-like,cnc,CrebA,CrebB,Dif"
txt="Atf3,BtbVII,Cf2,CG44247,CHES-1-like,crc,CrebB,Dif,dysf,E5,Ets98B,FoxP,Hnf4,Hr3,Hsf,jim,Mef2,Myc,onecut,pan,peb,Rbp6,sr,Trf2,trx,Xbp1"
txt="Atf3,bin,BtbVII,Cf2,CHES-1-like,cnc,CrebA,CrebB,Dif,foxo,ham,Hand,hb,Hnf4,jim,kay"
txt="Atf3,cnc,crc,CrebB,E2f1,foxo,GATAe,ham,Hnf4,kay,Mitf,pnt,Sox14,SREBP,vri,Xbp1"
txt="bon,RunxA"
txt="Doc3"
txt="BEAF-32,Hsf,Jra,bs,BtbVII,cg,CG4730,crc,Dp,Gsc,CG17829,maf-S,Max,slou,Trf2"
#txt="bs,Cf2,CG16779,CG7101,CHES-1-like,CrebB,E2f2,foxo,jim"
txt="Dfd,Dref,Ets98B,fkh,foxo,H15,Hnf4,Mef2,onecut,peb,Rbp6,slp1,sr,Trf2"
txt="Cf2,CHES-1-like,crc,E2f1,hng1,jim,Jra,Kdm4B,Mef2,Mondo,Myc,nej,ovo,Pdp1,srp,svp,Trf2,usp,vri"
txt="cnc,Trf2"
txt="toe"
txt="Hey,bcd,srl"
txt="dimm,Lim3"
txt="CG16779,CG4328,Lim3,onecut,otp,Psi,sr"
txt="dimm"
strings = txt.split(',')
tfgenenamedict1 = tfgenenamedict.copy()
# dict.fromkeys(['X', 'Y', 'Z'], 0)
print("===========rgenenamedict1=============")
# 초기화..
tfgenenamedict1 = tfgenenamedict1.fromkeys(tfgenenamedict1, 0)
print(tfgenenamedict1)
print(strings)
for i in strings:
    print(i)
    tfgenenamedict1[i] = 1
    alist.append(i)

print(tfgenenamedict1)
list1 = list(tfgenenamedict1.values())

print("========tf============")
#print(list1)
nparry = np.array([list1])
print("========tf nparry[0]============")
#print(nparry[0])

ry=rgenenamemodel.predict(nparry, verbose=0)


#X_train, X_test, y_train, y_test = train_test_split(arrtraining, labeltname, test_size=0.2,
#                                                    random_state=1)
print("========classification_report(rgene_test,predictedrgene========")
print(classification_report(rgene_test,predictedrgene))
print("=========cls_rpt===========")
#print(cls_rpt[['recall']])
print(rgene_test)
print(predictedrgene)

'''
X_train, X_test, y_train, rgene_test 
y_pred_keras = rgenenamemodel.predict(X_test).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred_keras)
test_predicted = rgenenamemodel.predict(X_test).ravel()
FPR, TPR, threshold = roc_curve(rgene_test, test_predicted)
plt.plot(FPR, TPR)
plt.plot([0,1],[0,1],'--')
plt.ylabel("True positive rate")
plt.xlabel("False positive rate")
'''

#ry=rgenenamemodel.predict(nparry)
print("==rgenenamemodel.predict===")
print(ry)

#y_predict = model.predict(x_test)
#print(y_predict) #[[0.67222077]]
#print(y_predict[0]) #[0.67222077]
#print(y_predict[0][0]) #0.67222077

#m= max(ry[0])
#print("=========m========")
#print(m)
#for i,x in enumerate(ry[0]):
#    print("index:{}, value :{}".format(i,x))





#최종 3개를 찾는것
n = 3
y_preds = np.argsort(ry, axis=1)[:,-n:]
print(y_preds)
print(rgenenamegetdict[int(y_preds[0][0])])
print(rgenenamegetdict[int(y_preds[0][1])])
print(rgenenamegetdict[int(y_preds[0][2])])
#최종 3개를 찾는것

print("tf=>",txt)
predictedrgene = ry.argmax(axis=-1)
print(predictedrgene)
print(rgenenamegetdict[int(predictedrgene)])


ty=tnamemodel.predict(nparry, verbose=0)

print("===========nparry==============")
print(nparry)
print("======ty=====")
print(ty)
predictedrgene = ty.argmax(axis=-1)
print("========predictedrgene==========")
print(predictedrgene)
print(tnamegetdict[int(predictedrgene)])




# initiate RMSprop optimizer  bad
# opt = keras.optimizers.rmsprop(lr=0.0001, decay=1e-6)
# model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])


'''
        [patterns.append(line.strip().split("\t")) for line in fp]
        #[patterns.append(line.strip().split(",")) for line in fp]
    text = f.read()
    words = text[:10000].split()
    stop_words = set(stopwords.words("english"))
    words_filter = [w for w in words if not w in stop_words]
              precision    recall  f1-score   support
           0       0.25      0.50      0.33         2
           1       1.00      0.88      0.93         8
           2       1.00      1.00      1.00         6
           3       0.67      1.00      0.80         8
           4       1.00      0.60      0.75        10
           5       1.00      1.00      1.00         2
           6       0.29      1.00      0.44         2
           7       1.00      1.00      1.00         3
           8       1.00      0.75      0.86         4
           9       0.00      0.00      0.00         1
          10       0.75      0.75      0.75         4
          11       0.80      0.80      0.80         5
          12       1.00      1.00      1.00         3
          13       0.83      1.00      0.91         5
          14       1.00      0.80      0.89         5
          15       1.00      0.67      0.80         3
          16       0.80      0.80      0.80         5
          17       1.00      1.00      1.00         4
          18       0.00      0.00      0.00         5
          19       0.44      1.00      0.62         4
          21       1.00      1.00      1.00         1
          22       1.00      1.00      1.00         3
          23       1.00      1.00      1.00         3
          24       1.00      0.60      0.75         5
          25       1.00      1.00      1.00         3
          26       0.62      0.83      0.71         6
          27       1.00      0.33      0.50         3
          28       1.00      0.90      0.95        10
          29       0.00      0.00      0.00         1
          30       0.50      0.50      0.50         2
          31       1.00      0.67      0.80         3
          32       1.00      0.50      0.67         2
          33       0.50      1.00      0.67         1
          34       0.83      1.00      0.91         5
          35       0.75      1.00      0.86         3
          36       1.00      1.00      1.00         5
          37       0.50      1.00      0.67         1
          40       0.67      0.67      0.67         3
          42       0.80      0.50      0.62         8
          43       1.00      1.00      1.00         4
          44       1.00      1.00      1.00         5
          45       0.50      0.62      0.56         8
          46       1.00      0.71      0.83         7
          47       1.00      1.00      1.00         3
          48       1.00      1.00      1.00         4
          49       0.00      0.00      0.00         1
          51       1.00      0.60      0.75         5
          52       1.00      0.50      0.67         2
          53       0.25      0.25      0.25         4
          54       0.83      0.62      0.71         8
          55       0.83      1.00      0.91         5
          56       0.57      0.80      0.67         5
          57       1.00      0.33      0.50         3
          58       1.00      0.67      0.80         3
          59       0.60      1.00      0.75         3
          60       1.00      1.00      1.00         3
          61       1.00      1.00      1.00         3
          62       0.33      0.50      0.40         2
          63       1.00      1.00      1.00         3
          64       1.00      0.50      0.67         4
          65       0.50      0.67      0.57         3
          66       0.67      1.00      0.80         4
          67       1.00      1.00      1.00         1
          68       1.00      1.00      1.00         2
          69       0.00      0.00      0.00         0
          70       1.00      1.00      1.00         3
          71       1.00      0.67      0.80         3
          72       1.00      0.40      0.57         5
          73       1.00      1.00      1.00         2
          74       0.71      1.00      0.83         5
          75       1.00      0.80      0.89         5
          77       1.00      1.00      1.00         5
          78       0.55      1.00      0.71         6
          79       1.00      1.00      1.00         3
          80       1.00      1.00      1.00         2
          81       0.78      0.88      0.82         8
          82       0.80      0.67      0.73         6

    accuracy                           0.79       305
   macro avg       0.79      0.77      0.75       305
weighted avg       0.83      0.79      0.79       305

'''
