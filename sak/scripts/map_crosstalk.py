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

print 'reading, %s'%opts.cal

#Array data -- why didn't grabbing the antpos from the cal file work? Got foobar.
aa = aipy.cal.get_aa(opts.cal, np.array([0.15]))
antpos = {
        '49': {'top_x':    0.0,  'top_y':    0.0,  'top_z':    0.0},
        '41': {'top_x': 3000.7,  'top_y':    6.5,  'top_z':    0.6},
        '47': {'top_x': 5999.1,  'top_y':    5.5,  'top_z':    6.9},
        '19': {'top_x': 8996.1,  'top_y':    8.1,  'top_z':    6.2},
        '29': {'top_x':11996.6,  'top_y':   12.2,  'top_z':    6.0},
        '28': {'top_x':14996.1,  'top_y':    5.2,  'top_z':   11.9},
        '34': {'top_x':17999.9,  'top_y':   10.9,  'top_z':    8.7},
        '51': {'top_x':20995.6,  'top_y':   15.7,  'top_z':   11.8},

        '10': {'top_x':   -1.8,  'top_y': -392.6,  'top_z':    8.8},
         '3': {'top_x': 2998.2,  'top_y': -394.1,  'top_z':    4.8},
        '25': {'top_x': 5998.9,  'top_y': -396.6,  'top_z':   13.5},
        '48': {'top_x': 8995.7,  'top_y': -393.4,  'top_z':   11.2},
        '24': {'top_x':12001.2,  'top_y': -390.7,  'top_z':   11.9},
        '55': {'top_x':14999.5,  'top_y': -392.1,  'top_z':   16.4},
        '27': {'top_x':17999.2,  'top_y': -386.9,  'top_z':   13.9},
        '57': {'top_x':20996.9,  'top_y': -384.1,  'top_z':   16.5},

         '9': {'top_x':   -1.1,  'top_y': -794.2,  'top_z':   11.8},
        '58': {'top_x': 2999.5,  'top_y': -793.5,  'top_z':    8.1},
         '1': {'top_x': 6000.9,  'top_y': -793.5,  'top_z':   16.5},
         '4': {'top_x': 8998.2,  'top_y': -795.0,  'top_z':   15.1},
        '17': {'top_x':12000.1,  'top_y': -790.6,  'top_z':   15.6},
        '13': {'top_x':14998.6,  'top_y': -793.4,  'top_z':   21.4},
        '56': {'top_x':17995.1,  'top_y': -789.7,  'top_z':   16.6},
        '59': {'top_x':20997.9,  'top_y': -785.8,  'top_z':   17.9},

        '22': {'top_x':   -3.9,  'top_y':-1193.5,  'top_z':   18.2},
        '61': {'top_x': 2999.3,  'top_y':-1191.8,  'top_z':   13.3},
        '35': {'top_x': 5999.3,  'top_y':-1194.8,  'top_z':   18.6},
        '18': {'top_x': 9000.4,  'top_y':-1193.4,  'top_z':   19.1},
         '5': {'top_x':11999.5,  'top_y':-1191.4,  'top_z':   17.8},
        '32': {'top_x':15000.4,  'top_y':-1193.2,  'top_z':   24.3},
        '30': {'top_x':18000.7,  'top_y':-1185.5,  'top_z':   20.6},
        '23': {'top_x':20998.6,  'top_y':-1185.9,  'top_z':   22.2},

        '20': {'top_x':   -2.0,  'top_y':-1591.7,  'top_z':   20.0},
        '63': {'top_x': 2999.3,  'top_y':-1593.8,  'top_z':   19.2},
        '42': {'top_x': 6000.5,  'top_y':-1592.1,  'top_z':   20.9},
        '37': {'top_x': 8997.8,  'top_y':-1589.1,  'top_z':   21.4},
        '40': {'top_x':11999.7,  'top_y':-1589.0,  'top_z':   20.7},
        '14': {'top_x':14998.5,  'top_y':-1590.5,  'top_z':   27.9},
        '54': {'top_x':17997.0,  'top_y':-1587.9,  'top_z':   22.8},
        '50': {'top_x':21002.3,  'top_y':-1587.3,  'top_z':   27.6},

        '43': {'top_x':   -1.8,  'top_y':-1994.9,  'top_z':   22.7},
         '2': {'top_x': 2999.9,  'top_y':-1994.0,  'top_z':   25.2},
        '33': {'top_x': 6000.7,  'top_y':-1991.0,  'top_z':   25.7},
         '6': {'top_x': 8999.8,  'top_y':-1990.2,  'top_z':   26.9},
        '52': {'top_x':11999.6,  'top_y':-1987.1,  'top_z':   23.9},
         '7': {'top_x':14996.1,  'top_y':-1986.6,  'top_z':   33.6},
        '12': {'top_x':17997.8,  'top_y':-1986.7,  'top_z':   28.6},
        '38': {'top_x':20998.2,  'top_y':-1983.3,  'top_z':   27.9},

        '53': {'top_x':   -1.8,  'top_y':-2393.3,  'top_z':   27.9},
        '21': {'top_x': 3000.7,  'top_y':-2392.3,  'top_z':   28.4},
        '15': {'top_x': 5998.7,  'top_y':-2391.6,  'top_z':   31.6},
        '16': {'top_x': 8999.9,  'top_y':-2389.6,  'top_z':   30.6},
        '62': {'top_x':11999.1,  'top_y':-2386.4,  'top_z':   28.0},
        '44': {'top_x':14996.5,  'top_y':-2387.0,  'top_z':   32.7},
         '0': {'top_x':17996.9,  'top_y':-2386.2,  'top_z':   31.4},
        '26': {'top_x':21000.5,  'top_y':-2386.2,  'top_z':   31.0},

        '31': {'top_x':   -0.6,  'top_y':-2796.2,  'top_z':   32.5},
        '45': {'top_x': 3001.4,  'top_y':-2793.6,  'top_z':   30.7},
         '8': {'top_x': 6001.6,  'top_y':-2790.5,  'top_z':   34.2},
        '11': {'top_x': 8998.1,  'top_y':-2790.8,  'top_z':   35.3},
        '36': {'top_x':12003.3,  'top_y':-2786.6,  'top_z':   32.0},
        '60': {'top_x':14997.5,  'top_y':-2786.5,  'top_z':   35.1},
        '39': {'top_x':18001.8,  'top_y':-2785.0,  'top_z':   33.2},
        '46': {'top_x':21001.2,  'top_y':-2784.8,  'top_z':   33.3},
}

#Recenter origin to (0,0) in bottom left
for i in range(64):
	antpos[str(i)]['top_x'] = antpos[str(i)]['top_x']-5.5
	antpos[str(i)]['top_y'] = antpos[str(i)]['top_y']+2795.6

#If wanted, plot the antenna grid
if opts.plotgrid:
	for i in range(64):
		vec = antpos[str(i)]
		pylab.plot(vec['top_x'],vec['top_y'],'ks')
		pylab.text(vec['top_x']+10.,vec['top_y'],str(i))
	pylab.show()

###
#Now we start focusing on a single antenna, and its relationships with the others
###
antorg = [antpos[opts.ant]['top_x'],antpos[opts.ant]['top_y']]

POS,DMAX,R,TH = [],[],[],[]


for i in range(64):
	if i==int(opts.ant): continue # no autos
	if i==37: continue #dammit zaki
	if i in badantarr: continue
	
	a = int(opts.ant)
	if a<i: bl='%i_%i'%(a,i)
	else: bl= '%i_%i'%(i,a)
	
	#Find the maximum crosstalk power in delay space
	W = P['counts'][bl]/np.float(np.max(P['counts'][bl]))
	weights = aipy.dsp.gen_window(P['counts'][bl].shape[0],'blackman-harris')
	_d,info = aipy.deconv.clean(np.fft.ifft(P[bl]*weights), np.fft.ifft(W),tol=1e-3)
	#delays = np.fft.fftfreq(len(W),d=np.diff(P['freqs'])[0])
	dmax_ij = np.amax(np.fft.fftshift(np.abs(_d+info['res']))) #<---
	
	other = [antpos[str(i)]['top_x'],antpos[str(i)]['top_y']]
	delta_x = antorg[0]-other[0]
	delta_y = antorg[1]-other[1]
	
	r = np.sqrt( np.power(delta_x,2) + np.power(delta_y,2) )
	
	#Still need to figure out why this happens sometimes. I'm avoiding autos! (1st line)
	try: theta = np.arctan(delta_y/delta_x)
	except ZeroDivisionError: continue
	
	if np.floor(np.abs(np.degrees(theta))) in range(88,92): continue
	
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
#Origin
axC.plot(antorg[0],antorg[1],'bs')
axC.text(antorg[0],antorg[1]+50.,opts.ant)
#Plot
nPOS = np.array(POS)
c = axC.scatter(nPOS[:,0], nPOS[:,1], s=DMAX)
pylab.show()
pylab.close()

if opts.plotrad:
	RD = np.array([R,DMAX])
	pylab.plot(RD[0],RD[1],'o')
	pylab.xlabel('Radial Distance')
	pylab.ylabel('Maximum xtalk in delay space')
	pylab.suptitle('Antenna %s'%opts.ant)
	pylab.show()	
	
	
	