#! /usr/bin/env python
import numpy as np, pylab as plt, capo
import sys

from keras.models import Sequential
from keras.layers import Convolution2D, MaxPooling2D
from keras.layers import Activation, Dropout, Flatten, Dense, Reshape

K = 8

model = Sequential()
model.add(Convolution2D(64, 3, 3, input_shape=(2*K, 2*K, 1), dim_ordering='tf', border_mode='same'))
model.add(Activation('relu'))
model.add(Flatten())
model.add(Dense((2*K)**2))
model.add(Activation('sigmoid'))
model.add(Reshape((2*K,2*K,1)))
model.add(Convolution2D(32, 3, 3, dim_ordering='tf', border_mode='same'))
model.add(Activation('relu'))
model.add(Flatten())
for cnt in xrange(5):
    model.add(Dense((5-cnt)*(2*K)**2))
    model.add(Activation('sigmoid'))
model.add(Reshape((2*K,2*K)))

model.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['accuracy'])

#model.load_weights('omni_xrfi_wgts_v000.hdf5')
#model.load_weights('omni_xrfi_wgts_v001.hdf5')
model.load_weights('omni_xrfi_wgts_v004.hdf5')

pol = 'xx'
X,Y = [],[]
for f in sys.argv[1:]:
    print 'Reading', f
    meta,g,v,xtalk = capo.omni.from_npz(f)
    w = meta['chisq']
    w_sm = np.empty_like(w)
    sig = np.empty_like(w)
    for i in xrange(w.shape[0]):
        for j in xrange(w.shape[1]):
            i0,j0 = max(0,i-K), max(0,j-K)
            i1,j1 = min(w.shape[0],i+K), min(w.shape[1],j+K)
            w_sm[i,j] = np.median(w[i0:i1,j0:j1])
    w_rs = w - w_sm
    w_sq = np.abs(w_rs)**2
    for i in xrange(w.shape[0]):
        for j in xrange(w.shape[1]):
            i0,j0 = max(0,i-K), max(0,j-K)
            i1,j1 = min(w.shape[0],i+K), min(w.shape[1],j+K)
            sig[i,j] = np.sqrt(np.median(w_sq[i0:i1,j0:j1]))
    f1 = w_rs / sig
    f1 = np.ma.array(f1, mask=np.where(f1 > 6, 1, 0))
    f1.mask |= np.isnan(f1)
    for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
        prev = 0
        x,y = np.where(f1.mask)
        while x.size != prev:
            prev = x.size
            xp, yp = (x+dx).clip(0,f1.shape[0]-1), (y+dy).clip(0,f1.shape[1]-1)
            i = np.where(f1[xp,yp] > 2)[0]
            f1.mask[xp[i],yp[i]] = 1
            x,y = np.where(f1.mask)
    f1ch = np.average(f1.mask, axis=0); f1ch.shape = (1,-1)
    f1.mask = np.logical_or(f1.mask, np.where(f1 > 6*(1-f1ch), 1, 0))
    f1t = np.average(f1.mask, axis=1)
    ts = np.where(f1t > 2*np.median(f1t))
    f1.mask[ts] = 1
    f1f_sum = np.sum(f1.filled(0), axis=0)
    f1f_wgt = np.sum(np.logical_not(f1.mask), axis=0)
    f1f = f1f_sum / f1f_wgt.clip(1,np.Inf)
    fs = np.where(f1f > 2)
    f1.mask[:,fs] = 1
    #capo.plot.waterfall(np.where(flg,0,w), drng=4)
    ##plt.subplot(131); capo.plot.waterfall(w, mx=0, drng=4)
    ##plt.subplot(132); capo.plot.waterfall(w_rs, mx=0, drng=4)
    ##plt.subplot(143); capo.plot.waterfall(sig, mx=0, drng=4)
    #capo.plot.waterfall(f1, mx=2, drng=2); plt.colorbar(); plt.show()
    #capo.plot.waterfall(f1.filled(10), mode='real', mx=10, drng=8); plt.show()
    #plt.show()
    #import IPython; IPython.embed()
    mask = f1.mask

    import IPython; IPython.embed()
    if False:
        for i in xrange(K,w.shape[0]-K,K/2):
            for j in xrange(K,w.shape[1]-K,K/2):
                i0,j0 = max(0,i-K), max(0,j-K)
                i1,j1 = min(w.shape[0],i+K), min(w.shape[1],j+K)
                #if np.random.uniform() > .1+ .9*mask[i0:i1,j0:j1].sum()/float(2*K)**2: continue
                X.append(np.random.uniform()*w[i0:i1,j0:j1])
                Y.append(mask[i0:i1,j0:j1])
X,Y = np.array(X), np.array(Y)
X.shape += (1,)
print 'Done'
import IPython; IPython.embed()
model.fit(X, Y, nb_epoch=15, batch_size=1000)

flg2_sum, flg2_wgt = np.zeros_like(w), np.zeros_like(w)
for i in xrange(K,w.shape[0]-K-1):
    print i
    for j in xrange(K,w.shape[1]-K-1,4):
        i0,j0 = i-K,j-K
        i1,j1 = i+K,j+K
        Xij = w[i0:i1,j0:j1]; Xij.shape = (1,2*K,2*K,1)
        Yij = model.predict(Xij)
        flg2_wgt[i0:i1,j0:j1] += 1
        flg2_sum[i0:i1,j0:j1] += Yij[0]
import IPython; IPython.embed()
flg2 = flg2_sum / flg2_wgt.clip(1,np.Inf)
wf1 = np.ma.array(w, mask=flg)
wf2 = np.ma.array(w, mask=np.where(flg2>.1,1,0))
plt.subplot(211); capo.plot.waterfall(wf1, mx=0, drng=4)
plt.subplot(212); capo.plot.waterfall(wf2, mx=0, drng=4)
plt.show()
wf1_rs = np.ma.array(w_rs/sig, mask=flg)
wf2_rs = np.ma.array(w_rs/sig, mask=np.where(flg2>.01,1,0))
plt.subplot(211); capo.plot.waterfall(wf1_rs, mx=2, drng=2)
plt.subplot(212); capo.plot.waterfall(wf2_rs, mx=2, drng=2)
plt.show()
g1 = np.ma.array(g['x'][88], mask=flg)
g2 = np.ma.array(g['x'][88], mask=np.where(flg2>.1,1,0))
plt.subplot(211); capo.plot.waterfall(g1, mx=0, drng=4)
plt.subplot(212); capo.plot.waterfall(g2, mx=0, drng=4)
plt.show()
