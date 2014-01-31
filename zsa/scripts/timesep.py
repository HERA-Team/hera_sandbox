#! /usr/bin/env python
import numpy as n
import optparse, sys, os

o = optparse.OptionParser()
o.add_option('--short', action='store_true',
            help='print short output. i.e. gaps or no gaps')
o.add_option('--file', action='store', default='timesep_outfile_spot.txt',
            help='filename to save dirs in')
o.add_option('--badfile', action='store', default='missing_files.txt',
             help='filename to save missing filenames')

opts,args = o.parse_args(sys.argv[1:])
if os.path.exists(opts.file):
    outfile = open(opts.file,'a') #write output of good files.
else: 
    outfile = open(opts.file, 'w')

if os.path.exists(opts.badfile):
    boutfile = open(opts.badfile,'a') #write output of bad files.
else: 
    boutfile = open(opts.badfile, 'w')

#if empty dir
if '*' in args[0]:
    boutfile.write(args[0].split('/')[0] + ' is an empty dir or does not exist.\n')
    print 'No dirs. Cleanly exiting'
    exit()
    

#script run on directories seperately.
#get first file and first time.
fin = args[0]
fsplit = fin.split('.') 
t0 = float(fsplit[1]+'.'+fsplit[2])
c = 0 
#time gap supposed tp be this. This is hard coded for ~10minute files.
sposed2 = 0.00696
#start loop on second file
for file in args[1:]:
    fsplit = file.split('.') 
    time = float(fsplit[1]+'.'+fsplit[2])
    diff = n.round(time-t0,decimals=5)
    if n.abs(n.round(diff - sposed2, decimals=5)) > 0.00001:
        if not opts.short:
            print 'gap between %s and %s is %f = %f files'%( fin.split('/')[-1],file.split('/')[-1], diff, diff/sposed2)
            for f in range(1,int(n.ceil(diff/sposed2)))[::-1]:
                boutfile.write('%s.%s.%d.%s\n'%(fsplit[0],fsplit[1],
                               int(fsplit[2])-(f*sposed2*1e5),fsplit[-1]))
        c+=1
    fin = file
    t0 = time
     
if c == 0: 
    print '%s has No Gaps! Ship It!'%file.split('/')[0]
    outfile.write(file.split('/')[0]+'\n')

if c>0: print 'Gaps. %s is a Bad Directory. %d gaps '%(file.split('/')[0], c)


