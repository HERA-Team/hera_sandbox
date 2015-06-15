#! /usr/bin/env python

import gspread
import os, sys
import optparse

o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('--email', dest='email', type='string',
            help='Google Spreadsheet Login Email (ex: name@gmail.com)')
o.add_option('--pw', dest='pw', type='string',  
            help='Google Spreadsheet Login Password')
o.add_option('--wkst', dest='wkst', type='string',
            help='Name of Google Spreadsheet Worksheet')
opts, args = o.parse_args(sys.argv[1:])

gc = gspread.login(opts.email,opts.pw)
wkst = gc.open(opts.wkst).sheet1

### HARD-CODED FOR 2014 PSA128 DATA ###
good_days_epoch1 = wkst.range('C2:C55')
good_days_epoch2 = wkst.range('D56:D116')
good_days_epoch3 = wkst.range('E117:E318')

file1 = open('good_days_epoch1.txt','w')
file2 = open('good_days_epoch2.txt','w')
file3 = open('good_days_epoch3.txt','w')

for c in range(len(good_days_epoch1)):
    val = good_days_epoch1[c].value
    file1.write(val)
    file1.write('\n')
file1.close()   
 
for c in range(len(good_days_epoch2)):
    val = good_days_epoch2[c].value
    file2.write(val)
    file2.write('\n')
file2.close()   
 
for c in range(len(good_days_epoch3)):
    val = good_days_epoch3[c].value
    file3.write(val)
    file3.write('\n')
file3.close()   

 
