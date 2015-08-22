#basic_aa.py
"""
This is a cal file for a basic amp antenna array arranged in a grid.
"""

import aipy as a, numpy as n,glob,ephem

prms = {
	'loc':('37:52:9.0','122:15:29.0',242), #What are the units for elevation? 242m
	'antpos':{
		# This is the large uv smear, small beam (sigma=0.25rad)
		0:[30./(4*n.pi),30./(4*n.pi),0.0],
		
		# These are for the small uv smears, large beams (sigma=1rad)
		1:[10./(2*n.pi),0.0,0.0],
		2:[10./(n.pi),0.0,0.0],
		3:[30./(2*n.pi),0.0,0.0],
		
		4:[0.0,10./(2*n.pi),0.0],
		5:[10./(2*n.pi),10./(2*n.pi),0.0],
		6:[10./(n.pi),10./(2*n.pi),0.0],
		7:[30./(2*n.pi),10./(2*n.pi),0.0],

		8:[0.0,10./(n.pi),0.0],
		9:[10./(2*n.pi),10./(n.pi),0.0],
		10:[10./(n.pi),10./(n.pi),0.0],
		11:[30./(2*n.pi),10./(n.pi),0.0],

		12:[0.0,30./(2*n.pi),0.0],
		13:[10./(2*n.pi),30./(2*n.pi),0.0],
		14:[10./(n.pi),30./(2*n.pi),0.0],
		15:[30./(2*n.pi),30./(2*n.pi),0.0],

		
	}
}

def get_aa(freqs):
	'''Return the AntennaArray to be used for the simulation.'''
	location = prms['loc']
	antennas = []
	nants = len(prms['antpos'])
	for ii in range(nants):
		pos = prms['antpos'][ii]
		beam = a.amp.Beam(freqs, xwidth=n.pi/8, ywidth=n.pi/8)
		antennas.append(a.amp.Antenna(pos[0],pos[1],pos[2],beam,phsoff=[0.,0.], 
			bp_r=n.array([1]), bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0)))
			#I'm just leaving all the kwargs as defaults since I'm not entirely sure what they are
	aa = a.amp.AntennaArray(location,antennas)
	return aa

def get_baselines_hybrid():
    na = len(aa.ants) # number of antennas
    baselines = n.zeros([16,3])
    for ii in n.arange(16):
        baselines[ii] = prms['antpos'][ii]
    return baselines
