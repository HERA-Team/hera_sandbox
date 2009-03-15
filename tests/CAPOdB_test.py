#!/usr/bin/env python
#
#  db_test.py
#  
#
#  Created by Danny Jacobs on 3/12/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import CAPOdB as cdb, sys, aipy as a,numpy as n

lat, lon = '45:00', '90:00'
freqs = a.loc.get_freqs(0.1,0.1/256,256)
beam = a.fit.BeamFlat(freqs)
aa = cdb.AntennaArray((lat,lon), [],beam=beam)
prms = aa.get_params()

aa[0].set_params(prms[0])
print "Testing initialization from paper_test1.db: ",
if prms[0]['x']==8.98:
   print "Passed"
else: print "Failed"; sys.exit()

print "Testing save/reload: ",
oldamp = aa[0].amp

aa[0].amp += 1
aa[0].save(aa.c,aa.aa_db_params)
aa[0].prm_update(aa.c,aa.aa_db_params)
prms2 = aa.get_params()
#print prms[0],prms2[0],aa[0].amp-oldamp-1
if n.abs(aa[0].amp-oldamp-1)<10**-10:
   print "Passed"
else: print "Failed"; sys.exit()
