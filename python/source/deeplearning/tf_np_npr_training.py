import logging
import os
import sys
import numpy as np
import pandas as pd
import pickle


import keras
from keras.datasets import cifar10
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
#from keras.optimizers import SGD
from keras.constraints import maxnorm


import tensorflow as tf


batch_size = 28
num_classes = 1001
epochs = 1000
data_augmentation = True
num_predictions = 1001
save_dir = os.path.join(os.getcwd(), 'saved_models')
imageSize=28



print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))  #D:\data\slb\project\python
basedir = os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))) #D:\data\slb\project\python
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
sys.path.append(os.path.dirname(os.path.dirname(__file__))) #load_info

import sysinfo
import platform
import pymysql


tfpath="D:/project/python/bioinfo/data/graph/tf"
modelpath = "D:/project/python/bioinfo/model"

filename="r_fca_biohub_testis_10x.loom"
filename="s_fca_biohub_antenna_10x.loom"
filename="s_fca_biohub_head_10x.loom"
filename="Chemosensory.txt"

npfilename="np_tf_tissue_training.txt"
npnprfilename="np_npr_tf_tissue_training.txt"
nprfilename="npr_tf_tissue_training.txt"
savetraining="save_np_npr_tf_tissue_training.pkl"

modeltname = "tname_training_tf"+str(epochs)+'.h5'
modelrgenename = "rgenename_training_tf"+str(epochs)+'.h5'

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




filedir= os.path.join(tfpath,npfilename)

dtnp = pd.read_csv(os.path.join(tfpath,npfilename), delimiter='\t')
print(len(dtnp))

dtnpr = pd.read_csv(os.path.join(tfpath,nprfilename), delimiter='\t')
print(len(dtnpr))

dtnpnpr = pd.read_csv(os.path.join(tfpath,npnprfilename), delimiter='\t')
print(len(dtnpnpr))

#dtfm['np_npr'] = dtfm['npgenename'].map(str) + '_' + dtfm['nprgenename'].map(str)
#nan 삭제..
#nprtdata = dtfm.dropna()
#column delete
#nprtdata =nprtdata.drop(['nprttype'],axis=1)


dtnpnpr_np = dtnpnpr.copy()
dtnpnpr_np= dtnpnpr_np.drop(['rgenename'],axis=1)
dtnpnpr_np.rename(columns = {'pgenename':'rgenename'},inplace=True)

print(dtnpnpr_np)

dtnpnpr_npr = dtnpnpr.copy()
dtnpnpr_npr= dtnpnpr_npr.drop(['pgenename'],axis=1)

print(dtnpnpr_npr)

dtnprtrain=pd.concat([dtnp,dtnpr,dtnpnpr_np,dtnpnpr_npr])
dtnprtrain.reset_index(inplace=True)

tnamedict ={}
npuniqtname= np.unique(dtnprtrain.tname)
print("======agenename[0]=========",npuniqtname[0])
for i in range(0,len(np.unique(dtnprtrain.tname))):
    tnamedict[npuniqtname[i]] = i

tfgenenamedict ={}
npuniqtfgenename= np.unique(dtnprtrain.tgenename)
print("======agenename[0]=========",npuniqtfgenename[0])
for i in range(0,len(np.unique(dtnprtrain.tgenename))):
    tfgenenamedict[npuniqtfgenename[i]] = i

rgenenamedict ={}
npuniqrgenename= np.unique(dtnprtrain.rgenename)
print("======agenename[0]=========",npuniqrgenename[0])
for i in range(0,len(np.unique(dtnprtrain.rgenename))):
    rgenenamedict[npuniqrgenename[i]] = i

print(tnamedict)
print(tfgenenamedict)
print(rgenenamedict)


print(dtnprtrain.tgenename)

tfgenenamedict ={}
alist =[]
for text in dtnprtrain.tgenename:
    strings = text.split(',')
    for i in strings:
        alist.append(i)
#print(alist)


npuniqtfgenename = np.unique(alist)
print(npuniqtfgenename)
for i in range(0,len(npuniqtfgenename)):
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

#print(dtnprtrain.rgenename[1])
#print(dtnprtrain.tname[0])


training=[]
icount=0

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


with open(os.path.join(tfpath,savetraining), 'rb') as save_training_file:
    # Step 3
    training = pickle.load(save_training_file)

    # After config_dictionary is read from file
    print(training)

print(tnamedict)
print(tfgenenamedict)
print(rgenenamedict)
#my_dict2 = dict((y,x) for x,y in my_dict.iteritems())

arrnp = np.array(training)

print(arrnp)
print(arrnp.shape)
print(len(arrnp))  #row
print(len(arrnp[0]))  #column

arrtraining = arrnp[:,0:450]
print(len(arrtraining))
print(len(arrtraining[0]))

labelrgenename =  arrnp[:,450]
labeltname =  arrnp[:,451]

print(labelrgenename.shape)
print(labeltname.shape)

#x_train, y_train , labelMap =load_plateimage.getDataLabel(input)
#x_test, y_test,labelMap =load_plateimage.getDataLabel(inputt)


# Convert class vectors to binary class matrices.
#y_train = keras.utils.to_categorical(y_train, num_classes)
#y_test = keras.utils.to_categorical(y_test, num_classes)




print("===============num_classes============")
print(num_classes)


print(arrtraining.shape[1:])

label = labelrgenename
#label = labeltname

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(arrtraining, label ,test_size=0.2,
                                                    random_state=1)
print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

print(X_train) #train value
print(X_test)
print(y_train) #train label
print(y_test)

num_classes = len(rgenenamedict)
#num_classes = len(tnamedict)


def create_model():
    model = keras.Sequential([
        keras.layers.Dense(64, activation='relu'),
        keras.layers.Dropout(0.25),
        keras.layers.Dense(64, activation='relu'),
        keras.layers.Dropout(0.25),
        keras.layers.Dense(num_classes, activation='softmax')
    ])



    model.compile(optimizer='adam',
                loss='sparse_categorical_crossentropy',
                metrics=['accuracy'])

    return model

model = create_model()

# Let's train the model using sgd  sgd is very good in now
lrate=0.001
decay=lrate/epochs
#sgd=SGD(lr=lrate, momentum=0.9, decay=decay, nesterov=False)
#model.compile(loss='categorical_crossentropy', optimizer=keras.optimizers.Adam(learning_rate=lrate) , metrics=['accuracy'])


'''
model.fit( arrtraining , labelrgenename,
          batch_size=batch_size,
          epochs=epochs,
          shuffle=True)
'''
#label = labeltname
label = labelrgenename
#label = labeltname

model.fit(arrtraining, label,
          batch_size=batch_size,
          epochs=epochs,
          shuffle=True)


# Save model and weights

modelname = "rgenename_training_tf"+str(epochs)+'.h5'
#modelname = "tname_training_tf"+str(epochs)+'.h5'


if not os.path.isdir(modelpath):
    os.makedirs(modelpath)
model_path = os.path.join(modelpath, modelname)
model.save(model_path)
print('Saved trained model at %s ' % model_path)



# initiate RMSprop optimizer  bad
#opt = keras.optimizers.rmsprop(lr=0.0001, decay=1e-6)
#model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])



'''
        [patterns.append(line.strip().split("\t")) for line in fp]
        #[patterns.append(line.strip().split(",")) for line in fp]
    text = f.read()
    words = text[:10000].split()
    stop_words = set(stopwords.words("english"))
    words_filter = [w for w in words if not w in stop_words]
    
    reconstructed_model = keras.models.load_model("my_model")
    
'''
