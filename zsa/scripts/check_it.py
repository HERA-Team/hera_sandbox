#!/usr/bin/env python
import glob

EXTENSION = '.uvcRREcAzxCPB'

fin = open('files_used_2_1_14', 'w')
c = 0
for f in glob.glob('psa*'):
    caled = len(glob.glob('%s/*_cal.npz'%f)) > 0
    used_in_lst = len(glob.glob('%s/*%s'%(f,EXTENSION))) > 0
    if caled and used_in_lst:
        print f + ' is caled and used.'
        c += 1
        fin.write(f+'\n')
        continue
    #if caled:
    #    print f + ' is caled.'
    #    continue
    #if used_in_lst:
    #    print f + ' is used in lst.'
    #    continue
    #    
    
print c
