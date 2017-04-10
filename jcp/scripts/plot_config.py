#! /usr/bin/env python

import numpy as n, pylab as p, sys
import matplotlib as mpl

hera = n.loadtxt('HERA352_UTM.dat')
paper = n.loadtxt('psa128.dat')

herasize = 14.6
papersize = 4.

#=============================== VERSION 0 ====================================
fig = p.figure()
ax = fig.add_subplot(111,aspect='equal')

#hera
for x,y in zip(hera[:,0],hera[:,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='tomato',radius=herasize/2.,fill=True))
    p.plot(x,y,color='white')
for x,y in zip(hera[:,0],hera[:,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='black',radius=herasize/2.,fill=False))
    p.plot(x,y,color='white')

#hera19
hera19ants = [0,1,2,11,12,13,14,23,24,25,26,27,37,38,39,40,52,53,54]
for x,y in zip(hera[hera19ants,0],hera[hera19ants,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='black',radius=herasize/2.,fill=True))
    p.plot(x,y,color='white')
hera19ants = [0,1,2,11,12,13,14,23,24,25,26,27,37,38,39,40,52,53,54]
for x,y in zip(hera[hera19ants,0],hera[hera19ants,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='cyan',radius=herasize/2.,fill=False))
    p.plot(x,y,color='white')


#hera37
hera37ants = [3,15,28,36,41,42,51,55,56,67,68,69,70,71,84,85,86,87]
for x,y in zip(hera[hera37ants,0],hera[hera37ants,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='gray',radius=herasize/2.,fill=True))
    p.plot(x,y,color='white')
hera37ants = [3,15,28,36,41,42,51,55,56,67,68,69,70,71,84,85,86,87]
for x,y in zip(hera[hera37ants,0],hera[hera37ants,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='orange',radius=herasize/2.,fill=False))
    p.plot(x,y,color='white')



p.xlim(540800,541150)
p.ylim(6.6e6 + 1050, 6.6e6 + 1350)
p.show()

sys.exit()

#=============================== VERSION 1 ====================================
fig = p.figure()
ax = fig.add_subplot(111,aspect='equal')

#hera
for x,y in zip(hera[:,0],hera[:,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='gray',radius=herasize/2.,fill=False,alpha=.5))
    p.plot(x,y,color='white')

#hera19
hera19ants = [0,1,2,11,12,13,14,23,24,25,26,27,37,38,39,40,52,53,54]
for x,y in zip(hera[hera19ants,0],hera[hera19ants,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='red',radius=herasize/2.,fill=False))
    p.plot(x,y,color='white')

#paper128
for x,y in zip(paper[:,0],paper[:,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='blue',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperhex1
paperhex1ants = [5,6,7,16,17,18,19,28,29,30,31,32,42,43,44,45,57,58,59]
for x,y in zip(hera[paperhex1ants,0],hera[paperhex1ants,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='purple',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperhex2
paperhex2ants = [142,143,144,162,163,164,165,182,183,184,185,186,202,203,204,205,221,222,223]
for x,y in zip(hera[paperhex2ants,0],hera[paperhex2ants,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='purple',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperoutriggers
paperoutriggers = [10,34,63,155,215,148,170,180,187,229,242,254,289,175,296,305,319,323]
for x,y in zip(hera[paperoutriggers,0],hera[paperoutriggers,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='green',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperrotate
paperrotate = [95,99,100,104,109,113,114,118,123,127]
for x,y in zip(paper[paperrotate,0],paper[paperrotate,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='red',width=papersize,height=papersize,fill=True))
    p.plot(x,y,color='white')

#paperstay
paperstay = [14,16,18,20,21,23,25,27,28,30,32,34,35,37,39,41,42,44,46,48,60,61,62,63,64,65,66,67,68,69,70,71,93,97,102,106,107,111,116,120,121,125]
for x,y in zip(paper[paperstay,0],paper[paperstay,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='blue',width=papersize,height=papersize,fill=True))
    p.plot(x,y,color='white')

p.xlim(540800,541150)
p.ylim(6.6e6 + 1050, 6.6e6 + 1350)
p.show()


#=============================== VERSION 2 ====================================
fig = p.figure()
ax = fig.add_subplot(111,aspect='equal')

#hera
for x,y in zip(hera[:,0],hera[:,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='gray',radius=herasize/2.,fill=False,alpha=.5))
    p.plot(x,y,color='white')

#hera19
hera19ants = [0,1,2,11,12,13,14,23,24,25,26,27,37,38,39,40,52,53,54]
for x,y in zip(hera[hera19ants,0],hera[hera19ants,1]):
    ax.add_artist(mpl.patches.Circle((x,y),color='red',radius=herasize/2.,fill=False))
    p.plot(x,y,color='white')

#paper128
for x,y in zip(paper[:,0],paper[:,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='blue',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperhex1
paperhex1ants = [5,6,7,16,17,18,19,28,29,30,31,32,42,43,44,45,57,58,59]
for x,y in zip(hera[paperhex1ants,0],hera[paperhex1ants,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='purple',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperrotate
paperrotate = [95,99,100,104,109,113,114,118,123,127]
for x,y in zip(paper[paperrotate,0],paper[paperrotate,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='red',width=papersize,height=papersize,fill=True))
    p.plot(x,y,color='white')

#paperstay
paperstay = [14,16,18,20,21,23,25,27,28,30,32,34,35,37,39,41,42,44,46,48,60,61,62,63,64,65,66,67,68,69,70,71,93,97,102,106,107,111,116,120,121,125]
for x,y in zip(paper[paperstay,0],paper[paperstay,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='blue',width=papersize,height=papersize,fill=True))
    p.plot(x,y,color='white')

#paperhexoutriggers
paperhexcenter = 30
hexdiam = 5*herasize
angle = 60*(n.pi/180) #60 degrees in radians
xoffsets1 = [-hexdiam*n.sin(angle),0,hexdiam*n.sin(angle)]
yoffsets1 = [hexdiam*n.cos(angle),hexdiam,hexdiam*n.cos(angle)]
xoffsets2,yoffsets2 = [],[]
for xoffset1 in xoffsets1:
    for xoffset2 in xoffsets1:
        xoffsets2.append(xoffset1 + xoffset2)
for yoffset1 in yoffsets1:
    for yoffset2 in yoffsets1:
        yoffsets2.append(yoffset1 + yoffset2)
xoffsets = xoffsets1 + xoffsets2 + [xoffsets1[2] + hexdiam*n.sin(angle)]
yoffsets = yoffsets1 + yoffsets2 + [yoffsets1[2] - hexdiam*n.cos(angle)]
for x,y in zip(xoffsets,yoffsets):
    x += hera[paperhexcenter,0] - papersize/2
    y += hera[paperhexcenter,1] - papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='orange',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperother
paperother = [3,4,15,36,41,50,51,55,56]
for x,y in zip(hera[paperother,0],hera[paperother,1]):
    x -= papersize/2
    y -= papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='magenta',width=papersize,height=papersize,fill=False))
    p.plot(x,y,color='white')

#paperimg attempt1: circle with jitter
center = 224
radius = 35. #meters
nants = 18
n.random.seed(1)
_xs,_ys = [],[]
cnt = 0
while True:
    theta = n.random.uniform(0,2*n.pi)
    xjitter = n.random.normal(scale=papersize/2.)
    yjitter = n.random.normal(scale=papersize/2.)
    x = radius * n.cos(theta) + xjitter
    y = radius * n.sin(theta) + yjitter
    #print theta,x,y
    #print n.array(_xs) - x, n.array(_ys - y) 
    if (n.abs(n.array(_xs) - x) < papersize).any() == True and (n.abs(n.array(_ys) - y) < papersize).any() == True: continue
    else: _xs.append(x),_ys.append(y)
    cnt += 1
    print cnt
    if cnt == nants: break    
#plotter
for x,y in zip(_xs,_ys):
    x += hera[center,0] - papersize/2
    y += hera[center,1] - papersize/2
    ax.add_artist(mpl.patches.Rectangle((x,y),color='green',width=papersize,height=papersize,fill=False))


p.xlim(540800,541150)
p.ylim(6.6e6 + 1050, 6.6e6 + 1350)
p.show()

