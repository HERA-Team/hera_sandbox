import sys,os,os.path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


prefix = 'v1_'
hexNums = [3, 6]
sep0 = 5
rots = [False,True]

#prefix = 'v2_'
#hexNums = [6, 3]
#sep0 = 5
#rots = [False,True]

#prefix = 'v3_'
#hexNums = [1, 24]
#sep0 = 8
#rots = [False,True]

prefix = 'v4_'
hexNums = [1, 30]
sep0 = 8
rots = [False,True]


#prefix = 'small_'
#hexNums = [3, 4]
#sep0 = 8
#rots = [False,True]


#prefix = 'v5_'
#hexNums = [2, 3, 3]
#sep0 = 3
#rots = [False,True,False]
#
#prefix = 'v6_'
#hexNums = [2, 2, 2, 2]
#sep0 = 3
#rots = [True,False,True,False]

SplitCore = False
SplitCoreOutriggers = 4
CutDownTo = None
SpreadOut128 = False
precisionFactor = 1000 

seps = [sep0]
for n in range(1,len(hexNums)): seps = [seps[0] * 2*(hexNums[n]-1)*3**.5] + seps

#right = sep0*np.asarray([1,0,0])
#upRight = sep0*np.asarray([.5,3**.5/2,0])
#upLeft = sep0*np.asarray([-.5,3**.5/2,0])

positions = [[0,0,0]]
for sep,hexNum,rot in zip(seps,hexNums,rots):
    newPos = []
    for pos in positions:
        for row in range(hexNum-1,-(hexNum)+SplitCore,-1):
            for col in range(0,2*hexNum-abs(row)-1):
                if rot:
                    xPos = pos[0] + ((-(2*hexNum-abs(row))+2)/2.0 + col)*sep;
                    yPos = pos[1] + row*sep*3**.5/2;
                else:
                    yPos = pos[1] + ((-(2*hexNum-abs(row))+2)/2.0 + col)*sep;
                    xPos = pos[0] + row*sep*3**.5/2;
                if prefix == 'v4_':
                    if (xPos**2+yPos**2 < 175**2): newPos.append([xPos, yPos, 0])
                else: newPos.append([xPos, yPos, 0])
    positions = newPos


    

if True:
    baselineDict = {}
    print len(positions)
    for ant1 in range(len(positions)):
        for ant2 in range(ant1+1,len(positions)):
            baseline = tuple([int(np.round(precisionFactor*(positions[ant1][i]-positions[ant2][i]))) for i in range(3)])
            if baselineDict.has_key(baseline): baselineDict[baseline] += 1
            else: baselineDict[baseline] = 1


    uniqueBaselines = np.asarray(baselineDict.keys())/(1.0*precisionFactor)
    redundancy = np.asarray(baselineDict.values())
    print "With", len(positions), "antennas there are", len(uniqueBaselines), "unique baselines."
    uniqueBaselines = np.append(uniqueBaselines, -uniqueBaselines, axis=0)
    redundancy = np.append(redundancy, redundancy, axis=0)
    uniqueBaselines = uniqueBaselines[np.argsort(redundancy)]
    redundancy = redundancy[np.argsort(redundancy)]

#%%
    plt.figure(1); plt.clf()
    plt.scatter(np.asarray(positions)[:,0],np.asarray(positions)[:,1],c=[.5,.5,.5], linewidths=.5, s=20)
    plt.xlabel('Position (m)'); plt.ylabel('Position (m)')
    plt.axis('square')
    plt.xlim(-200,200); plt.ylim(-200,200)
    plt.savefig(prefix+'heir.pdf')

    plt.figure(2); plt.clf()
#    plt.scatter(uniqueBaselines[:,0]/1.0/1.5, uniqueBaselines[:,1]/1.0/1.5,s=7,c=(np.minimum(redundancy,20)),edgecolors='face')
    plt.scatter(uniqueBaselines[:,0]/1.0/1.5, uniqueBaselines[:,1]/1.0/1.5,s=7,c=np.log10(redundancy),cmap='inferno',edgecolors='face')
    plt.xlabel('Baseline (m)'); plt.ylabel('Baseline (m)')
#    plt.title('log10(redundancy)')
    plt.colorbar(label='$\log_{10}($Simultaneous Redundancy$)$')
    plt.axis('square')
    plt.xlim(-300,300); plt.ylim(-300,300)
    plt.savefig(prefix+'base.png', dpi=300)
#%%
else:
    plt.figure(1, figsize=(8,8)); plt.clf()
    plt.scatter(np.asarray(positions)[:,0],np.asarray(positions)[:,1])
    plt.xlabel('Position (m)'); plt.ylabel('Position (m)')
    plt.axis('square')

if False:                
    import omnical.info as oi
    import omnical.arrayinfo as oai
    redundantInfo = oi.RedundantInfo()
    positions = np.asarray(positions)
    reds = oai.compute_reds(positions)
    redundantInfo.init_from_reds(reds, positions)
    
    AtransA = redundantInfo.At.dot(redundantInfo.At.T).toarray()
    BtransB = redundantInfo.Bt.dot(redundantInfo.Bt.T).toarray()    
    gainVariances = np.diag(np.linalg.pinv(AtransA)[0:len(positions),0:len(positions)])
    print "Gain and phase modes that can't be Omnicaled:"
    print [len(XtransX) - np.linalg.matrix_rank(XtransX, tol=1e-10) for XtransX in [AtransA, BtransB]]
#%%
    plt.figure(3); plt.clf()
    hexScatter = plt.scatter(positions[:,0], positions[:,1], c=np.sqrt(gainVariances), cmap='inferno', linewidths=.5)
    plt.axis('square')
    #plt.xlim(-200,200); plt.ylim(-200,200)
    plt.colorbar(label='Antenna Relative Gain Errors')
    
    plt.xlabel('Position (m)'); plt.ylabel('Position (m)')
    plt.xlim(-200,200); plt.ylim(-200,200)
    plt.savefig(prefix+'cal.pdf')

plt.show()