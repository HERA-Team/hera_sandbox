#basic_aa.py
"""
This is a cal file for a basic antenna array.
"""

import aipy as a, numpy as n,glob,ephem

prms = {
	'loc':('37:52:9.0','122:15:29.0',242), #What are the units for elevation? 242m
	'antpos':{
		0:[0.0,0.0,0.0],
		1:[1.0,100.0,0.0],
		2:[2.0,10.0,5.0],
		3:[10.0,50.0,10.0]
	}
}

def get_aa(freqs):
	'''Return the AntennaArray to be used for the simulation.'''
	location = prms['loc']
	antennas = []
	nants = len(prms['antpos'])
	for ii in range(nants):
		pos = prms['antpos'][ii]
		beam = a.phs.Beam(freqs)
		antennas.append(a.phs.Antenna(pos[0],pos[1],pos[2],beam,phsoff=[0.,0.])) #What is phsoff?
	aa = a.phs.AntennaArray(location,antennas)
	return aa

