#! /usr/bin/env python
import WatchAlert
import numpy as np
import matplotlib.pyplot as plt

w = WatchAlert.Watch()
a = WatchAlert.Alert()
b = WatchAlert.Bandwidth()

skip = True
if not skip:

    watchFiles = ['watch_140903.txt','watch_141013.txt','watch_141107.txt','watch_141120.txt','watch_150115.txt','watch_150128.txt','watch_150721.txt','watch_151015.txt']

    for wf in watchFiles:
        w.readFile(wf)
    w.writeFile()

    a.readFile()
    vvv = np.ones(len(a.alerts))*370

    plt.plot(w.times,w.ping,'.',label='ping')
    plt.plot(a.alerts,vvv,'s',label='alerts')
    plt.ylabel('ms')
    plt.legend()


b.readbwFile()
bw = b.realtime/b.filesize*1024.*1024.*8

plt.figure('Bandwidth')
plt.plot(b.times, bw, 'k.')
#plt.plot(b.times, bw)
plt.ylabel('Mbps')
plt.show()