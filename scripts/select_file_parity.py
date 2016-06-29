import sys,optparse,re,os
"""
select files with even or odd jd
"""
o = optparse.OptionParser()
o.set_usage('nightly_avg.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--odd',action='store_true',
    help="Select odd days. (default is even)")
opts,args = o.parse_args(sys.argv[1:])

def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])
for filename in args:
    jd_int = int(file2jd(os.path.basename(filename)))
    if jd_int%2 and opts.odd:
        print filename
    elif jd_int%2==0 and not opts.odd:
        print filename
