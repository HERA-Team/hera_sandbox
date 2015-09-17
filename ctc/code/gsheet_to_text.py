#! /usr/bin/env python

import gspread
import json
from oauth2client.client import SignedJwtAssertionCredentials
import os, sys
import optparse

o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('--OAuth', dest='OAuth', type='string',
            help='Path to .json Google OAuth code.')
o.add_option('--wkst', dest='wkst', type='string',
            help='Name of Google Spreadsheet Worksheet')
opts, args = o.parse_args(sys.argv[1:])

json_key = json.load(open(opts.OAuth))
scope = ['https://spreadsheets.google.com/feeds']
credentials = SignedJwtAssertionCredentials(json_key['client_email'],json_key['private_key'],scope)
gc = gspread.authorize(credentials)
wkst = gc.open(opts.wkst).sheet1

### HARD-CODED FOR 2014 PSA128 DATA ###
good_days_epoch1 = wkst.range('C2:C55')
good_days_epoch2 = wkst.range('C56:C116')
good_days_epoch3 = wkst.range('C117:C246')

file1 = open('good_days_epoch1.txt','w')
file2 = open('good_days_epoch2.txt','w')
file3 = open('good_days_epoch3.txt','w')

for c in range(len(good_days_epoch1)):
    val = good_days_epoch1[c].value
    file1.write(val)
    file1.write('\n')
print file1,'written'
file1.close()   
 
for c in range(len(good_days_epoch2)):
    val = good_days_epoch2[c].value
    file2.write(val)
    file2.write('\n')
print file2,'written'
file2.close()   
 
for c in range(len(good_days_epoch3)):
    val = good_days_epoch3[c].value
    file3.write(val)
    file3.write('\n')
print file3,'written'
file3.close()   

 
