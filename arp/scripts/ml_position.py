#! /usr/bin/env python

import aipy as ap, numpy as np, pylab as plt, capo
import sys

from keras.models import Sequential
from keras.layers import Convolution1D, MaxPooling1D
from keras.layers import Activation, Dropout, Flatten, Dense
from keras.preprocessing.image import ImageDataGenerator, array_to_img, img_to_array, load_img

NCHAN = 1024

model = Sequential()
model.add(Activation('relu', input_shape=(NCHAN,1)))
model.add(MaxPooling1D(pool_length=4, border_mode='same'))

model.add(Convolution1D(32, 16, border_mode='same'))
model.add(Activation('relu'))
model.add(MaxPooling1D(pool_length=4, border_mode='same'))

model.add(Convolution1D(32, 16, border_mode='same'))
model.add(Activation('relu'))
model.add(MaxPooling1D(pool_length=4, border_mode='same'))

model.add(Flatten())  # this converts our 3D feature maps to 1D feature vectors
model.add(Dense(64))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(1))
model.add(Activation('sigmoid'))

model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])
try: model.load_weights('sep01.hdf5')
except(IOError): pass

# get some data
info,data,flgs = capo.miriad.read_files(sys.argv[1:], antstr='cross', polstr='xx,yy')
POL = data.values()[0].keys()[0]
aa = ap.cal.get_aa('hsa7458_v000_HH', np.array([.15]))
info = capo.omni.aa_to_info(aa)
reds = info.get_reds()
sep10 = [r for r in reds if (112,72) in r][0]
bad_ant = [81]
X,Y = [],[]
for bl in data:
    if bad_ant in bl: continue
    Xi = np.abs(data[bl][POL]) # don't use phase information, only fringe interference
    if bl in sep10 or bl[::-1] in sep10:
        Yi = np.ones(Xi.shape[0])
    else:
        Yi = np.zeros(Xi.shape[0])
    #X.append(Xi[:5]); Y.append(Yi[:5])
    X.append(Xi); Y.append(Yi)
X,Y = np.concatenate(X, axis=0), np.concatenate(Y, axis=0)
X.shape += (1,)
model.fit(X, Y, nb_epoch=50, batch_size=1000)
import IPython; IPython.embed()
    
