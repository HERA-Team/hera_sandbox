import pylab,sys

filestr = sys.argv[1:][0]
print filestr

F = open(filestr,'r')
T,RMS = [],[]
label,badlist,goodlist = [],[],[]

rmsbound = 0.2

c = 0
print 'BAD:'
for line in F: 
	L = line.split()
	#if L[0] in badlist_02: continue
	#if L[0] in badlist_0075: continue
	T.append(float(L[0][4:17]))
	RMS.append(float(L[1]))
	label.append(L[0][0:17])
	
	if float(L[1]) >= rmsbound: 
		print float(L[0][4:17]), float(L[1])
		badlist.append(L[0])
	else: goodlist.append(L[0])
	c+=1
pylab.plot(T,RMS,'-',lw=2)
pylab.axhline(rmsbound,linestyle='--',color='k')

pylab.margins(0.2)

pylab.show()
pylab.close()

print 'GOOD:',len(goodlist),'/',c,'rms-bound of %f'%rmsbound
