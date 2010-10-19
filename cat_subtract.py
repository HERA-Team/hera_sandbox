#!/usr/bin/env python
#
#  cat_subtract.py
#  
#
#  Created by Danny Jacobs on 3/18/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem,logging,os,warnings,vo.table as vot
from astrogrid import ConeSearch
from cStringIO import StringIO
"""
Subtract catalog A from catalog B.  ie output a new catalog C such that none
of the entries listed in B remain.

Usage
cat_subtract.py -s all --cat=catA,catB
"""



o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('--sep',dest='sep',type='float',default=0.25,
    help="Maximum acceptable seperation in degrees [0.25]")
o.add_option('-v',dest='verb',action='store_true',
    help='Print more')    
o.add_option('-f',dest='freq',default=0.15,type='float',
    help='Frequency at which to generate catalog')
o.add_option('--voenhance',dest='vo',action='store_true',
    help="Add alternate names, specfind id if available.  Requires vo desktop to be running.")
o.add_option('--uri',default="ivo://CDS.VizieR/VIII/85",
    help="URI for vo data source that supports cone search. default=ivo://CDS.VizieR/VIII/85 [specfindv2]")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")
opts, args = o.parse_args(sys.argv[1:])


def find_src_in_cat(src,cat,sep):
    """
    Given a src (src = aipy.amp.RadioFixedObject or lower) 
    and a catalog (cat = aipy.amp.RadioCatalog or lower), find the closest
    source in cat that falls within sep [deg].
    Return: closest source object and seperation (ephem.Angle).
    If no source found within sep, return None,closest_sep
    """
    last_sep=n.pi
    closest_src=None
    for nk,test_src in cat.iteritems():
        cursep = ephem.separation(src, test_src)
        if cursep <= last_sep:
            last_sep = cursep
            closest_src = test_src
    if last_sep<sep * a.img.deg2rad: return closest_src,last_sep
    else: return None,last_sep

def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)

def format_src(src):
    out = src.src_name
    out += '\t' + str(src._ra) 
    out += '\t' + str(src._dec)
    out += '\t' + str(src._jys)
    return out

def unit(votarray,item):
    units = {'mJys':1/1000.,
            'mJy':1/1000.}
    for i,t in enumerate(votarray.array.dtype.names):
        if t.startswith(item):uname=votarray.fields[i].unit
    try: return units[uname]
    except(KeyError): 1.

def find_src_name(votarray):
    name_prefs = ('3C','VLSS','NVSS','TXS')
    for cat_name in name_prefs:
        try:
            for name in votarray.array['Name']:
                if name.startswith(cat_name): return name
        except(KeyError): return "no names"
        except AttributeError:
            for name in votarray.array['Name']:
                if name.startswith(cat_name): return name        
        
    return votarray.array['Name'][0]    
def find_seq_of_brightest_src(vodata):
    return list(set(vodata[n.where(vodata['S_nu_']==n.max(vodata['S_nu_']))]['Seq']))[0]

#setup logging
if opts.vverb: logging.disable(logging.DEBUG)
elif opts.verb: logging.disable(logging.INFO)
else: 
    logging.disable(logging.ERROR)
    warnings.simplefilter('ignore',Warning)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)


#warnings.simplefilter('ignore',n.RankWarning)

catA_name = opts.cat.split(',')[0]
catB_name = opts.cat.split(',')[1]

if not opts.cal is None:
    catA = a.cal.get_catalog(opts.cal,catalogs=[catA_name])
    catB = a.cal.get_catalog(opts.cal,catalogs=[catB_name])
else:
    catA = a.src.get_catalog(catalogs=[catA_name])
    catB = a.src.get_catalog(catalogs=[catB_name])
    
update_pos(catA)
catA.update_jys(opts.freq)
update_pos(catB)
catB.update_jys(opts.freq)
print "#"+' '.join(sys.argv)
if opts.vo:
    print "#Name \t Ra \t Dec \t S_nu_ \t e_S_nu_ \t Altnames \t Seq"
else:
     print "#Name \t Ra \t Dec \t S_nu_ \t e_S_nu_"
    

for name,src in catA.iteritems():
    fsrc,sep = find_src_in_cat(src,catB,opts.sep)
    if fsrc is not None: continue
    try: err=src.e_S_nu
    except AttributeError: err = 0
#    print "%s \t %s \t"%(src.ra,src.dec),
    #perform the cone search and make a spectrum numpy array
    if opts.vo:
        log.debug("Connecting to %s"%(opts.uri))
        conesrch = ConeSearch(opts.uri)
        log.debug("Performing conesearch")
        votable = conesrch.execute(src.ra/ephem.degree,src.dec/ephem.degree,opts.sep)
        log.debug("Converting votable string to file object")
        votable_f = StringIO(votable)
        log.debug("Parsing table from votable to numpy.")
        votable_np = vot.parse_single_table(votable_f,pedantic=False)
        vodata = votable_np.array
        if len(vodata)<2: 
            print "%s \t %s \t %s \t %07.3f \t %07.3f"%\
                (name,src.ra,src.dec,src.jys,err)
            log.warning("no vo sources found within %7.3f of %s"%(opts.sep,name))
            continue
        #extract a single spectrum (hope there was only one hit, throw a warning at the end if not)
        seqs = list(set(vodata['Seq']))
        if len(seqs)>1: log.warning("%d vo sources found within %7.3f of %s"%(len(seqs),opts.sep,name))
        seq = find_seq_of_brightest_src(vodata)
        vodata = vodata[n.argwhere(vodata['Seq']==seq)]
        altnames = ','.join(vodata['Name'].squeeze())
        print "%s \t %s \t %s \t %07.3f \t %07.3f \t %s \t %d"%\
                (name,src.ra,src.dec,src.jys,err,altnames,int(seq))
    else:
        print "%s \t %s \t %s \t %07.3f \t %07.3f"%\
            (name,src.ra,src.dec,src.jys,err)        

        


