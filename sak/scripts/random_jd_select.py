import sys,os,optparse,re
import numpy as np
"""

"""
o = optparse.OptionParser()
o.add_option('--seed', type='int', default=10, help='Random seed. Choose one!')
o.add_option('--choice', action='store_true', help="Which bin's worth of filenames to be printed? Bin1 <=> True, Bin2 <=> False.")
o.add_option('--verbose', action='store_true', help='Toggle verbosity.')
opts,args = o.parse_args(sys.argv[1:])

def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])

#get integer JDs
ijds = []
for filename in args: ijds.append(int(file2jd(os.path.basename(filename))))
ijdlist = list(set(ijds))

#we are splitting the list in half
if not len(ijdlist)%2 == 0: S = len(ijdlist)+1
else: S = len(ijdlist)

#set random seed if provided
if not opts.seed is None: np.random.seed(opts.seed)

#split JDs
c1 = np.random.choice(ijdlist,size=S/2,replace=False)
c2 = np.array(filter(lambda jd: jd not in c1, ijdlist))

if opts.verbose:
    print " ".join(map(str,sorted(c1)))
    print " ".join(map(str,sorted(c2)))

for filename in args:
    ijd = int(file2jd(os.path.basename(filename)))
    if opts.choice:
        if ijd in c1:
            print filename
    else:
        if ijd in c2:
            print filename
