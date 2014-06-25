#basic_aa.py
"""
This is a cal file for a basic amp antenna array arranged in a circle.
"""

import aipy as a, numpy as n,glob,ephem

prms = {
	'loc':('37:52:9.0','122:15:29.0',242), #What are the units for elevation? 242m
	'antpos':{
		# this arrangement lowered the difference to 10 (from 100 for paper)
		0:[0.0,0.0,0.0],

		# radius 2ns
		1:[2.0,0.0,0.0],
		2:[0.0,2.0,0.0],
		3:[-2.0,0.0,0.0],
		4:[0.0,-2.0,0.0],

		5:[1.4,1.4,0.0,0.0],
		6:[-1.4,1.4,0.0,0.0],
		7:[-1.4,1.4,0.0,0.0],
		8:[1.4,-1.4,0.0,0.0],
		
		# radius 5ns -- gah... these are not circles, they're stars....
		# 1:[5.0,0.0,0.0],
		# 2:[0.0,5.0,0.0],
		# 3:[-5.0,0.0,0.0],
		# 4:[0.0,-5.0,0.0],

		# 5:[7.0,7.0,0.0,0.0],
		# 6:[-7.0,7.0,0.0,0.0],
		# 7:[-7.0,7.0,0.0,0.0],
		# 8:[7.0,-7.0,0.0,0.0],

		# radius 10 ns
		# 1:[10.0,0.0,0.0],
		# 2:[0.0,10.0,0.0],
		# 3:[-10.0,0.0,0.0],
		# 4:[0.0,-10.0,0.0],

		# 5:[14.0,14.0,0.0,0.0],
		# 6:[-14.0,14.0,0.0,0.0],
		# 7:[-14.0,14.0,0.0,0.0],
		# 8:[14.0,-14.0,0.0,0.0],

		
	}
}

def get_aa(freqs):
	'''Return the AntennaArray to be used for the simulation.'''
	location = prms['loc']
	antennas = []
	nants = len(prms['antpos'])
	for ii in range(nants):
		pos = prms['antpos'][ii]
		beam = a.amp.Beam(freqs, poly_azfreq=n.array([[.5]]))
		antennas.append(a.amp.Antenna(pos[0],pos[1],pos[2],beam,phsoff=[0.,0.], 
			bp_r=n.array([1]), bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0)))
			#I'm just leaving all the kwargs as defaults since I'm not entirely sure what they are
	aa = a.amp.AntennaArray(location,antennas)
	return aa

