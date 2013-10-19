#! /usr/bin/env python
import sys, os
import aipy as a

cat2ned = {
    'parkes':'pks',
    'culgoora':'Cul',
}

def ned_curl(src, cat):
    os.system('curl "http://ned.ipac.caltech.edu/cgi-bin/datasearch?objname=%s+%s&meas_type=bot&ebars_spec=ebars&label_spec=no&x_spec=freq&y_spec=Fnu_jy&xr=-1&of=ascii_bar&search_type=Photometry" > /tmp/%s.spec' % (cat, src, src))
    return open('/tmp/%s.spec' % src).read()

for srclist in sys.argv[1:]:
    srclist = [src for src in open(srclist).read().split()]
    print srclist
    for catalog in ['parkes','culgoora']:
        nedcat = cat2ned[catalog]
        cat = a.src.get_catalog(srclist, catalogs=[catalog])
        for src in cat:
            filename = '%s_spec.txt' % src
            data = ned_curl(src, nedcat).splitlines()
            print 'Writing', filename
            f = open(filename,'w')
            f.write('\n'.join(data[17:])) # skip 17 lines for header in NED
            f.close()
