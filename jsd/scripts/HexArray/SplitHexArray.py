import sys,os,os.path
import numpy as np
import cPickle as pickle
import matplotlib
import matplotlib.pyplot as plt
import pylab

Separation = 14.6
hexNum = 11
SplitCore = True
SplitCoreOutriggers = 0
CutDownTo = 242
SpreadOut128 = False
precisionFactor = 1000 

#Calculating Positions:
#Main Hex
positions = [];
for row in range(hexNum-1,-(hexNum)+SplitCore,-1):
#for row in range(hexNum-1,-(hexNum),-1):
    #for col in range(2*hexNum-abs(row)-1):
    for col in range(0,2*hexNum-abs(row)-1):
        xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*Separation;
        yPos = row*Separation*3**.5/2;
        positions.append([xPos, yPos, 0])


right = Separation*np.asarray([1,0,0])
up = Separation*np.asarray([0,1,0])
upRight = Separation*np.asarray([.5,3**.5/2,0])
upLeft = Separation*np.asarray([-.5,3**.5/2,0])

def cutSide(positions, th, saven=None, saveth=None):
    positions = np.asarray(positions)
    rotator = np.asarray([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0,0,1]])
    rotatedPositions = rotator.dot(positions.T).T
    maxy =  np.max(rotatedPositions[:,1])
    toCut = np.abs(rotatedPositions[:,1] - maxy) <.01
    if saven is not None:
        rotator = np.asarray([[np.cos(saveth), -np.sin(saveth), 0], [np.sin(saveth), np.cos(saveth), 0], [0,0,1]])
        rotatedPositions = rotator.dot(positions.T).T
        toCut = np.logical_and(toCut, rotatedPositions[:,1] < sorted(rotatedPositions[toCut][:,1])[-saven])
    return list(positions[np.logical_not(toCut),:])

#HERA 128 Spread Out
if len(positions)==320 and CutDownTo == 128 and SpreadOut128:
       
    coreHex = np.asarray(positions)
    coreHex = cutSide(coreHex, 0)
    coreHex = cutSide(coreHex, np.pi/3)
    coreHex = cutSide(coreHex, -np.pi/3, saven=1, saveth=np.pi)
    coreHex = cutSide(coreHex, 2*np.pi/3)
    coreHex = cutSide(coreHex, np.pi/3)  
    coreHex = cutSide(coreHex, 0)
    coreHex = cutSide(coreHex, 0)
    coreHex = cutSide(coreHex, np.pi/3)       
    
    for i in range(4):
        positions = cutSide(positions, 0)
        positions = cutSide(positions, np.pi/3)
        positions = cutSide(positions, -np.pi/3)
        positions = cutSide(positions, 2*np.pi/3)
        positions = cutSide(positions, np.pi/3)  
        positions = cutSide(positions, 0)
    positions = cutSide(positions, 0)
    positions = cutSide(positions, 0)
    positions = cutSide(positions, np.pi/3)
    positions = cutSide(positions, np.pi/3)
    positions = cutSide(positions, 2*np.pi/3)
    
 
    topHex = np.asarray(coreHex)
    for i in range(6):
        topHex = cutSide(topHex, np.pi/3)
        topHex = cutSide(topHex, 2*np.pi/3)
        topHex = cutSide(topHex, 3*np.pi/3)
        topHex = cutSide(topHex, 4*np.pi/3)
        topHex = cutSide(topHex, 2*np.pi/3)
        topHex = cutSide(topHex, 3*np.pi/3)        
    topHex = cutSide(topHex, np.pi/3)
    topHex = cutSide(topHex, 2*np.pi/3)
    topHex = cutSide(topHex, 2*np.pi/3)

    rightHex = np.asarray(coreHex)
    for i in range(5):
        rightHex = cutSide(rightHex, 5*np.pi/3)
        rightHex = cutSide(rightHex, 5*np.pi/3)
        rightHex = cutSide(rightHex, 4*np.pi/3)
        rightHex = cutSide(rightHex, 4*np.pi/3)
        rightHex = cutSide(rightHex, 3*np.pi/3)
        rightHex = cutSide(rightHex, 0*np.pi/3)        
    for i in range(2): rightHex = cutSide(rightHex, 0*np.pi/3)
    for i in range(3): rightHex = cutSide(rightHex, 4*np.pi/3)
    for i in range(2): rightHex = cutSide(rightHex, 5*np.pi/3)
  
        
    positions = np.vstack((positions, topHex, rightHex))

#HERA 128 Compact
if len(positions)==320 and CutDownTo == 128 and not SpreadOut128:
    for i in range(3):
        positions = cutSide(positions, 0)
        positions = cutSide(positions, np.pi/3)
        positions = cutSide(positions, -np.pi/3)
        positions = cutSide(positions, 2*np.pi/3)
        positions = cutSide(positions, np.pi/3)  
        positions = cutSide(positions, 0)
    positions = cutSide(positions, -np.pi/3)#, saven=1, saveth=0)
    positions = cutSide(positions, 2*np.pi/3)
    positions = cutSide(positions, 2*np.pi/3, saven=3, saveth=0)

#HERA 243
if len(positions)==320 and CutDownTo == 243:
    positions = cutSide(positions, 0)
    positions = cutSide(positions, np.pi/3)
    positions = cutSide(positions, -np.pi/3, saven=1, saveth=np.pi)
    positions = cutSide(positions, 2*np.pi/3)
    positions = cutSide(positions, np.pi/3)  
    positions = cutSide(positions, 0)
    positions = cutSide(positions, 0)
    positions = cutSide(positions, np.pi/3)
    
#HERA 37
if len(positions)==320 and CutDownTo == 19:
    for i in range(6):
        positions = cutSide(positions, 0)
        positions = cutSide(positions, np.pi/3)
        positions = cutSide(positions, -np.pi/3)
        positions = cutSide(positions, 2*np.pi/3)
        positions = cutSide(positions, np.pi/3)  
        positions = cutSide(positions, 0)
    positions = cutSide(positions, 2*np.pi/3)
    positions = cutSide(positions, 2*np.pi/3)
    positions = cutSide(positions, 0)
    positions = cutSide(positions, np.pi/3)
    positions = cutSide(positions, np.pi/3)
    CutDownTo = None

originalPositions = np.asarray(positions)

if len(positions)==320 and CutDownTo == 242:
    positions = cutSide(positions, 0)
    positions = cutSide(positions, np.pi/3)
    positions = cutSide(positions, -np.pi/3)#, saven=1, saveth=np.pi)
    positions = cutSide(positions, 2*np.pi/3)
    positions = cutSide(positions, np.pi/3)  
    positions = cutSide(positions, 0)
    positions = cutSide(positions, np.pi)
    positions = cutSide(positions, 4*np.pi/3)




if CutDownTo is not None:
    while len(positions) > CutDownTo:                                
        positions = cutSide(positions, 0)           
        if len(positions) <= CutDownTo: break
        positions = cutSide(positions, np.pi/3)            
        if len(positions) <= CutDownTo: break
        positions = cutSide(positions, -np.pi/3) 
        if len(positions) <= CutDownTo: break
        positions = cutSide(positions, 2*np.pi/3)            
        if len(positions) <= CutDownTo: break
        positions = cutSide(positions, np.pi/3)
        if len(positions) <= CutDownTo: break
        positions = cutSide(positions, 0)
       

#Split the core into 3 pieces
if SplitCore:
    newPos = []
    for i,pos in enumerate(positions):          
        theta = np.arctan2(pos[1],pos[0])
        if (pos[0]==0 and pos[1]==0):
            newPos.append(pos)
        elif (theta > -np.pi/3 and theta < np.pi/3):
            newPos.append(np.asarray(pos) + (upRight + upLeft)/3)                    
        elif (theta >= np.pi/3 and theta < np.pi):
            newPos.append(np.asarray(pos) +upLeft  - (upRight + upLeft)/3)
        else:
            newPos.append(pos)
    positions = newPos

nCore = len(positions)

if SplitCoreOutriggers:
    exteriorHexNum = SplitCoreOutriggers
    for row in range(exteriorHexNum-1,-(exteriorHexNum),-1):
        for col in range(2*exteriorHexNum-abs(row)-1):
            xPos = ((-(2*exteriorHexNum-abs(row))+2)/2.0 + col)*Separation*(hexNum-1)
            yPos = row*Separation*(hexNum-1)*3**.5/2
            theta = np.arctan2(yPos,xPos)       
            if ((xPos**2 + yPos**2)**.5 > Separation*(hexNum+1)):
                if (theta > 0 and theta <= 2*np.pi/3+.01):
                    positions.append(np.asarray([xPos, yPos, 0]) - 4*(upRight + upLeft)/3)
                elif (theta <= 0 and theta > -2*np.pi/3):
                    positions.append(np.asarray([xPos, yPos, 0])- 2*(upRight + upLeft)/3)
                else:
                    positions.append(np.asarray([xPos, yPos, 0]) - 3*(upRight + upLeft)/3)



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


plt.figure(1, figsize=(18,8)); plt.clf()
plt.subplot(121)
plt.scatter(np.asarray(positions)[:,0],np.asarray(positions)[:,1])
plt.xlabel('Position (m)'); plt.ylabel('Position (m)')
plt.axis('square')

plt.subplot(122)
plt.scatter(uniqueBaselines[:,0]/1.0/1.5, uniqueBaselines[:,1]/1.0/1.5,s=7,c=(np.minimum(redundancy,20)),cmap='inferno',edgecolors='face')
#plt.scatter(uniqueBaselines[:,0]/1.0/1.5, uniqueBaselines[:,1]/1.0/1.5,s=7,c=np.log10(redundancy),cmap='inferno',edgecolors='face')
plt.xlabel('Baseline (m)'); plt.ylabel('Baseline (m)')
plt.title('log10(redundancy)')
plt.xlim(-300,300); plt.ylim(-300,300)

#%%
#Saving results
scriptDirectory = os.path.dirname(os.path.abspath(__file__))
np.savetxt(scriptDirectory + "/antenna_positions.dat",np.asarray(positions))
plt.figure(2)
for n in [350,243,128,37]:
    array = np.loadtxt('./ProposedFinalConfigs/antenna_positions_'+str(n)+'.dat')
    plt.plot(array[:,0], array[:,1],'o')
plt.xlabel('Position (m)'); plt.ylabel('Position (m)')
plt.savefig('/Users/jsdillon/Desktop/HERA_Array_Stages.pdf')

# if __name__ == "__main__":
#     #from mpldatacursor import datacursor        
#     uniqueBaselines = np.asarray([uniqueBaseline[0] for uniqueBaseline in baselineDict.items()])/(1.0*precisionFactor)
#     redundancy = np.asarray([len(uniqueBaseline[1]) for uniqueBaseline in baselineDict.items()])
#     uniqueBaselines = np.append(uniqueBaselines, -uniqueBaselines, axis=0)
#     redundancy = np.append(redundancy, redundancy, axis=0)
#     plt.figure(1)
#     plt.clf()
#     #plt.plot(originalPositions[:,0]/Separation, originalPositions[:,1]/Separation,'r.')        
#     plt.scatter(np.asarray(positions)[:,0],np.asarray(positions)[:,1])
    
#     plt.axis('equal')
#     plt.figure(2)    
#     plt.clf()
# #        plt.scatter(uniqueBaselines[:,0]/1.0/Separation, uniqueBaselines[:,1]/1.0/Separation,c=np.minimum(redundancy,100000),s=40)
#     plt.scatter(uniqueBaselines[:,0]/1.0/1.5, uniqueBaselines[:,1]/1.0/1.5,s=40,c=np.log10(np.minimum(redundancy,40)))
#     plt.xlim(-350,350); plt.ylim(-350,350)
# #        plt.colorbar()
#     plt.title('Baseline Redundancy')
#     #datacursor(display='single',formatter="x={x:.4f}\ny={y:.4f}".format)
#     #plt.axis('equal')
#     return np.asarray(positions), np.asarray(positions[0:nCore])

# if __name__ == "__main__":
#     positions,corePos = HexArray(hexNum = 11, SplitCore = True, SplitCoreOutriggers = True, CutDownTo = None, SpreadOut128 = False)

if True:                
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
    plt.figure(3)
    plt.clf()
    hexScatter = plt.scatter(positions[:,0], positions[:,1], c=np.sqrt(gainVariances), s=100)
    plt.colorbar()
    plt.title('Antenna Relative Gain Errors')
#%%
if True:
    import omnical.info as oi
    import omnical.arrayinfo as oai
    corePos = positions[0:nCore,:] 
    coreRedundantInfo = oi.RedundantInfo()
    coreReds = oai.compute_reds(corePos)
    coreRedundantInfo.init_from_reds(coreReds, corePos)
    coreAtransA = coreRedundantInfo.At.dot(coreRedundantInfo.At.T).toarray()
    coreBtransB = coreRedundantInfo.Bt.dot(coreRedundantInfo.Bt.T).toarray()    
    coreGainVariances = np.diag(np.linalg.pinv(coreAtransA)[0:len(corePos),0:len(corePos)])
    print "Gain and phase modes that can't be Omnicaled:"
    print [len(XtransX) - np.linalg.matrix_rank(XtransX, tol=1e-10) for XtransX in [coreAtransA, coreBtransB]]


#     if True:
# #        plt.figure(101); plt.clf()        
# #        for array in [np.loadtxt('antenna_positions_330.dat'), np.loadtxt('antenna_positions_243_toward_330.dat')]:
# #            plt.plot(array[:,0], array[:,1],'o')
#         plt.figure(101); plt.clf()               
# #        for n in [350,243,128,37]:
#         for n in [350,243,128,37]:
#             array = np.loadtxt('antenna_positions_'+str(n)+'_alt.dat')
#             plt.plot(array[:,0], array[:,1],'o')
# #        plt.plot(positions[:,0], positions[:,1], 'ko')
#%%
if True:
    plt.figure(10); plt.clf()
    bls = {bl[1]: red for bl,red in zip(uniqueBaselines, redundancy) if bl[0] == 0 and bl[1] > 0}
    freqs = np.arange(50,250,200.0/256)
    temp = plt.scatter(bls.values(),bls.keys(), c=np.log10(bls.values()), cmap="inferno"); plt.clf()
    for bl,red in bls.items():
        if bl < 500: 
            plt.plot(bl*freqs * 1e6/3e8,freqs, c=plt.cm.get_cmap("inferno")(np.log10(red)/np.log10(np.max(bls.values()))))
    plt.colorbar(temp, label='$\log_{10}($Simultaneous Redundancy$)$')
    plt.xlabel('|u| (wavelengths)'); plt.ylabel('Frequency (MHz)')
    plt.savefig('HERA Hex u Coverage.pdf')

plt.show()

#%% For Proposal

if True:
    def plotArray(antennaPositions):
        fig = plt.gcf()
        xLim = 1.2 * np.max(antennaPositions)
        yLim = xLim
        plt.xlim(-xLim,xLim)
        plt.ylim(-yLim,yLim)
        for antennaPos in antennaPositions:
            fig.gca().add_artist(plt.Circle((antennaPos[0],antennaPos[1]),7,fc='0.3',ec='k'))
    
    def dishScatter(fig, ax, xPos, yPos, cVals, colormap, radii=7.0, cLim = None, blackOutlines = 0): 
        if blackOutlines: radii -= .5
        if not hasattr(radii,'len'): radii = np.ones(len(xPos))*radii    
        if cLim is None:     
            underlyingScatter = plt.scatter(xPos, yPos, c=cVals, s=0, cmap=colormap.name)
            cVals = (cVals - np.min(cVals))/(np.max(cVals) - np.min(cVals))
        else:
            underlyingScatter = plt.scatter(xPos, yPos, c=cVals, s=0, cmap=colormap.name, vmin = cLim[0], vmax = cLim[1])
            cVals = 1.0*(cVals - cLim[0])/(cLim[1]-cLim[0])
        for x,y,r,c in zip(xPos, yPos, radii, colormap(cVals)): 
            ax.add_patch(plt.Circle((x,y), r, fc=c, ec='k',lw=blackOutlines))
        return underlyingScatter

    
#%%
    positions = np.asarray(positions)
    import matplotlib.gridspec as gridspec
    fig = plt.figure(101, figsize=(14*1.3,5*1.3)); plt.clf()
    gs = gridspec.GridSpec(1,3, width_ratios=[1,1,1])
    plt.subplots_adjust(hspace = 0.0, wspace = .25, bottom = .24, top=.99, left = .042, right = .98)
    axes = [plt.subplot(gs[i]) for i in range(3)]

    hexScatter = dishScatter(fig, axes[0], positions[:,0], positions[:,1], np.ones(len(positions)), pylab.cm.gray, 7.0, cLim=[.6,2.0], blackOutlines = .25)    
    axes[0].set_ylim([-140,140])        
    axes[0].set_xlim([-145,135])         
    axes[0].set_xlabel('Position (m)',size=16)
    axes[0].set_ylabel('Position (m)',size=16, labelpad=-10)
    axes[0].set_aspect('equal', adjustable='box')
    

    #allBaselines = np.append(np.loadtxt('unique_baselines.dat'),-np.loadtxt('unique_baselines.dat'),axis=0)
    allRedundancies =redundancy# np.append(np.loadtxt('redundancy.dat'),np.loadtxt('redundancy.dat'),axis=0)
    allBaselines = uniqueBaselines
    scatterPlot = axes[1].scatter(allBaselines[:,0],allBaselines[:,1],c=np.log10(allRedundancies),edgecolors='none',s=22,vmin=0, vmax=2.5, cmap=pylab.cm.inferno, linewidths=.25)
    axes[1].set_ylim([-260,260])        
    axes[1].set_xlim([-260,260])        
    axes[1].set_ylabel('Baseline (m)',size=14,labelpad=-10)            
    axes[1].set_xlabel('Baseline (m)',size=14)            
    axes[1].set_aspect('equal', adjustable='box')
    
    
#        thisScatter = dishScatter(fig, axes[2], corePos[:,0], corePos[:,1], np.sqrt(coreGainVariances), pylab.cm.inferno, 7.0, blackOutlines = .5)
    corePos = positions[0:nCore,:] 
    thisScatter = dishScatter(fig, axes[2], corePos[:,0], corePos[:,1], np.sqrt(coreGainVariances), pylab.cm.inferno, 7.0, blackOutlines = .5, cLim=[0.0653, .0703])
    axes[2].set_xlim([-145,135])        
    axes[2].set_ylim([-140,140])    
    axes[2].set_xlabel('Position (m)',size=12)
    axes[2].set_ylabel('Position (m)',size=12,labelpad=-10)
    axes[2].set_aspect('equal', adjustable='box')

    cax = plt.axes([axes[1].get_position().x0, .1, axes[1].get_position().x0 + axes[1].get_position().width - axes[1].get_position().x0, .04])
    clr = plt.colorbar(scatterPlot, cax=cax,orientation='horizontal', ticks=np.arange(0,3,.5))      
    clr.set_label('$\log_{10}($Simultaneous Redundancy$)$', size=16)              


    cax2 = plt.axes([axes[2].get_position().x0, .1, axes[2].get_position().x0 + axes[2].get_position().width - axes[2].get_position().x0, .04])
    clr2 = plt.colorbar(thisScatter, cax=cax2, orientation='horizontal', ticks=[.066, .067, .068, .069, .07])                    
    clr2.set_label('Relative Gain Calibration Error', size=16) 
    
    for ax in axes:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(16)    
    
    plt.savefig('/Users/jsdillon/Desktop/HERA_Array_Config.pdf', format='pdf')