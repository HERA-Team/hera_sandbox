"""
The PAPER catalog here is derived by experimental methods in either facets or 
healpix maps.
"""

import aipy as a, numpy as n, os,logging,vo.table as vot

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('generic_catalog')

class RadioFixedBody(a.fit.RadioFixedBody):
    def __init__(self,*args,**kwargs):
        a.fit.RadioFixedBody.__init__(self, *args, **kwargs)
        if kwargs.has_key('e_S_nu'):self.e_S_nu = kwargs['e_S_nu']
        if kwargs.has_key('vo_record'): self.vo_record = kwargs['vo_record']

    
class GENERICCatalog(a.fit.SrcCatalog):
    def fromfile(self, filename):
        addsrcs = []
        data = n.genfromtxt(filename,names=True,dtype=None)
        altnames = {'S_nu':'jys','nu':'mfreq','Ra':'ra','Dec':'dec'}
        prms = {'jys':0,'Name':'','index':0,'mfreq':0.15,'e_S_nu':0,'ra':0,'dec':0}
        for rec in data:
            for prm in data.dtype.names:
                try:
                    if prms.has_key(prm): prms[prm]=rec[prm]
                    elif prms.has_key(altnames[prm]):
                        prms[altnames[prm]] = rec[prm]
                    else: 
                        log.error("column not recognized")
                except(KeyError):
                        log.warning("Skipping unknown parameter "+prm)
            addsrcs.append(RadioFixedBody(prms['ra'],prms['dec'],name=prms['Name'],
                jys=prms['jys'],index=prms['index'],mfreq=prms['mfreq'],e_S_nu=prms['e_S_nu']))
        self.add_srcs(addsrcs)
                    
#    def fromfile(self, filename):
#        f = open(filename)
#        addsrcs = []
#        for L in [L for L in f.readlines() if not L.startswith('#')]:
#            text = L.split('\t')
#            if len(text) <= 3: continue
##            try: int(text[0].strip()[0])
##            except(ValueError): continue
#            ra = text[1]
#            dec = text[2]
#            name = text[0].strip()
#            jys = float(text[3])
#            try: e_S_nu = float(text[4])
#            except(IndexError): e_S_nu = 0
#            addsrcs.append(RadioFixedBody(ra, dec, name=name,
#                jys=jys, index=0, mfreq=.150,e_S_nu=e_S_nu))
#        self.add_srcs(addsrcs)
    def fromvot(self,filename):
        addsrcs = []
        votable = vot.parse_single_table(filename,pedantic=False)
        for rec in votable.array:
            ra = rec['Ra']
            dec = rec['Dec']
            name = rec['Name']
            jys = rec['S_nu_']
            try: e_S_nu = rec['e_S_nu_']
            except(IndexError): e_S_nu = 0
            try: index = rec['n']
            except(IndexError):
                index = rec['index']
            addsrcs.append(RadioFixedBody(ra, dec, name=name,
                jys=jys, index=index, mfreq=.150,e_S_nu=e_S_nu,vo_record=rec))            
        self.add_srcs(addsrcs)

#GENERICFILE = os.path.dirname(__file__) + os.sep + 'generic.txt'
_genericcat = {}

def get_srcs(srcs=None, cutoff=None,catalogs=None):
    global _genericcat
    if catalogs is None: return []
    if type(catalogs) is list: catalog=catalogs[0]
    else: catalog = catalogs
    GENERICFILE = os.path.dirname(__file__) + os.sep+catalog+'.txt'
    GENERICVOT = os.path.dirname(__file__) + os.sep+catalog+'.vot'
    try: 
        open(GENERICVOT)
        GENERICFILE=GENERICVOT
        log.info("opening %s"%(GENERICFILE))
    except(IOError): 
        try: 
            open(catalog+'.vot')
            GENERICFILE=catalog+'.vot'
            log.info("opening %s"%(GENERICFILE))
        except(IOError): 
            try: 
                open(GENERICFILE)
                log.info("opening %s"%(GENERICFILE))
            except(IOError):
                try: 
                    open(catalog+'.txt')
                    GENERICFILE=catalog+'.txt'
                    log.info("opening %s"%(GENERICFILE))
                except(IOError):
                    log.warning("%s %s"%(GENERICFILE,"not found"))
                    return []
    if not _genericcat.has_key(GENERICFILE):
        log.info("reading file")
        _genericcat[GENERICFILE] = GENERICCatalog()
        if GENERICFILE.endswith('.txt'):
            _genericcat[GENERICFILE].fromfile(GENERICFILE)
        elif GENERICFILE.endswith('.vot'):
            _genericcat[GENERICFILE].fromvot(GENERICFILE)
        log.info("loaded %d sources from %s"%(len(_genericcat[GENERICFILE]),GENERICFILE))
    if srcs is None:
        if cutoff is None: srcs = _genericcat[GENERICFILE].keys()
        else:
            cut, fq = cutoff
            fq = n.array([fq])
            for s in _genericcat[GENERICFILE].keys(): _genericcat[GENERICFILE][s].update_jys(fq)
            srcs = [s for s in _genericcat[GENERICFILE].keys() if _genericcat[GENERICFILE][s].jys[0] > cut]

    srclist = []
    for s in srcs:
        try:
            try: srclist.append(_genericcat[GENERICFILE][s.src_name])
            except(AttributeError):
                srclist.append(_genericcat[GENERICFILE][s])
        except(KeyError): pass
    return srclist
