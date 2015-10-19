"""
This script plots the output of uv_avg_SK.py or nightly_avg_SK.py (SK=polarization options, see capo/sak/scripts, otherwise see capo/dcj/scripts). 

Written by: Saul A. Kohn, Daniel C. Jacobs
"""

import numpy as np, sys, optparse, ephem, os, pylab, pickle, aipy

o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('map_crosstalk.py [options] xtalk.pkl')
aipy.scripting.add_standard_options(o, cal=True)
o.add_option('-a','--ant',dest='ant',default='49',help='antenna to plot xtalk around.')
o.add_option('-b','--badant',dest='badant',default=None,help='bad antennae to remove, separated by commas. e.g. "-b 1,2,3"')
o.add_option('--plotgrid',dest='plotgrid',action='store_true',help='Plot the antenna grid positions and display before doing other stuff.')
o.add_option('--plotrad',dest='plotrad',action='store_true',help='Make a plot of delay peak vs radial distance from chosen antenna.')

opts, args = o.parse_args(sys.argv[1:])

#Parse bad antennae
badantarr = []
if opts.badant != None:
	nocomma = opts.badant.split(',')
	for entry in nocomma:
		badantarr.append(int(entry))


#TODO: Right now it can only take one pkl file. I should relax this restriction and have it loop through and save separate figs. 
print 'loading, %s'%args[0]

P = pickle.load(open(args[0]))
freqs = P['freqs']*1e3

#Array data
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']

Nants = len(antpos.keys())

#Recenter origin to (0,0) in bottom left
for i in range(Nants):
	if Nants==64: #PSA-64 grid
		antpos[i]['top_x'] = antpos[i]['top_x']-5.5
		antpos[i]['top_y'] = antpos[i]['top_y']+2795.6
	if Nants==128: #PSA-128 grid
		antpos[i]['top_x'] = antpos[i]['top_x']-540875+0.2
		antpos[i]['top_y'] = antpos[i]['top_y']-6601130-3.5-0.2
		
		
#If wanted, plot the antenna grid
if opts.plotgrid:
	for i in range(Nants):
		vec = antpos[i]
		pylab.plot(vec['top_x'],vec['top_y'],'ks')
		pylab.text(vec['top_x'],vec['top_y'],str(i))
	pylab.show()

###
#Now we start focusing on a single antenna, and its relationships with the others
###
antorg = [antpos[int(opts.ant)]['top_x'],antpos[int(opts.ant)]['top_y']]

POS,DMAX,R,TH = [],[],[],[]


for i in range(Nants):
	if i==int(opts.ant): continue # no autos
	if Nants==64 and i==37: continue #dammit zaki
	if i in badantarr: continue
	
	a = int(opts.ant)
	if a<i: bl='%i_%i'%(a,i)
	else: bl= '%i_%i'%(i,a)
	
	#Find the maximum crosstalk power in delay space
	if np.float(np.max(P['counts'][bl])) != 0.: W = P['counts'][bl]/np.float(np.max(P['counts'][bl]))
	else: 
		print 'no counts for BL=%s?'%bl
		continue
	weights = aipy.dsp.gen_window(P['counts'][bl].shape[0],'blackman-harris')
	_d,info = aipy.deconv.clean(np.fft.ifft(P[bl]*weights), np.fft.ifft(W),tol=1e-3)
	#delays = np.fft.fftfreq(len(W),d=np.diff(P['freqs'])[0])
	
	#put whatever stat you like in here. I've put the maximum point in delay space there.
	if Nants==64: dmax_ij = np.amax(np.fft.fftshift(np.abs(_d+info['res'])))
	elif Nants==128: dmax_ij = 1e5*np.amax(np.fft.fftshift(np.abs(_d+info['res'])))
	
	other = [antpos[i]['top_x'],antpos[i]['top_y']]
	delta_x = antorg[0]-other[0]
	delta_y = antorg[1]-other[1]
	
	r = np.sqrt( np.power(delta_x,2) + np.power(delta_y,2) )
	
	#Still need to figure out why this happens sometimes. I'm avoiding autos! (1st line)
	try: theta = np.arctan(delta_y/delta_x)
	except ZeroDivisionError: continue
	
	if np.floor(np.abs(np.degrees(theta))) in range(88,92): continue
	
	#Store calculations in lists
	POS.append(other)
	DMAX.append(dmax_ij)
	R.append(r)
	TH.append(theta)
	#axC.plot([antorg[0],other[0]],[antorg[1],other[1]],'b-') <-- "asterix" plot


#axP = pylab.subplot(111,polar=True) <-- why doesn't this work? angles aren't working properly...

#Plot the xtalk/grid scatter by default. Should I make this optional? Perhaps just the .show() optional and save the fig?

#Setup
axC = pylab.subplot(111)
axC.set_xlabel('E-W')
axC.set_ylabel('N-S')
axC.set_title('Antenna %s'%opts.ant)
#Origin
axC.plot(antorg[0],antorg[1],'rs')
axC.text(antorg[0],antorg[1]+50.,opts.ant)
#Plot
nPOS = np.array(POS)
c = axC.scatter(nPOS[:,0], nPOS[:,1], s=DMAX)
pylab.show()
pylab.close()

if opts.plotrad:
	RD = np.array([R,DMAX])
	
	par = np.polyfit(np.log10(R), np.log10(DMAX), 1, full=True)
	slope=par[0][0]
	intercept=par[0][1]
	print 'log(y) propto',slope,'log(R)'
	
	pylab.loglog(RD[0],RD[1],'o')
	
	if Nants==64:
		R30 = np.sqrt( pow(antpos[49]['top_x']-antpos[41]['top_x'],2) + pow(antpos[49]['top_y']-antpos[41]['top_y'],2) )
		R60 = np.sqrt( pow(antpos[49]['top_x']-antpos[47]['top_x'],2) + pow(antpos[49]['top_y']-antpos[47]['top_y'],2) )
		R90 = np.sqrt( pow(antpos[49]['top_x']-antpos[19]['top_x'],2) + pow(antpos[49]['top_y']-antpos[19]['top_y'],2) )
		
		for Rx in [R30,R60,R90]: pylab.axvline(Rx,color='k',ls=':')
		
	if Nants==128: 
		R15 = np.sqrt( pow(antpos[64]['top_x']-antpos[10]['top_x'],2) + pow(antpos[64]['top_y']-antpos[10]['top_y'],2) )
		R30 = np.sqrt( pow(antpos[64]['top_x']-antpos[49]['top_x'],2) + pow(antpos[64]['top_y']-antpos[49]['top_y'],2) )
		R45 = np.sqrt( pow(antpos[64]['top_x']-antpos[3]['top_x'],2) + pow(antpos[64]['top_y']-antpos[3]['top_y'],2) )
		R60 = np.sqrt( pow(antpos[64]['top_x']-antpos[41]['top_x'],2) + pow(antpos[64]['top_y']-antpos[41]['top_y'],2) )
		R75 = np.sqrt( pow(antpos[64]['top_x']-antpos[25]['top_x'],2) + pow(antpos[64]['top_y']-antpos[25]['top_y'],2) )
		R90 = np.sqrt( pow(antpos[64]['top_x']-antpos[19]['top_x'],2) + pow(antpos[64]['top_y']-antpos[19]['top_y'],2) )
		
		for Rx in [R15,R30,R45,R60,R75,R90]: pylab.axvline(Rx,color='k',ls=':')
	
	#xl = [np.amin(R), np.amax(R)]
	#yl = [slope*xx + intercept  for xx in xl]
	
	pylab.xlabel('Radial Distance')
	pylab.ylabel('Maximum xtalk in delay space')
	pylab.suptitle('Antenna %s'%opts.ant)
	pylab.show()	
	

	