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
            help='Name of Google Spreadsheet Worksheet.')
o.add_option('--pg', dest='pagenum', type='int',
            help='Page number of Spreadsheet.')
o.add_option('-s', dest='season', type='int',
            help='Season of 128 Data on Google Doc (1 or 2).')
o.add_option('-p',dest='pol',type='string',default='xx',
            help='Pol of data. Default is xx.')
opts, args = o.parse_args(sys.argv[1:])

json_key = json.load(open(opts.OAuth))
scope = ['https://spreadsheets.google.com/feeds']
credentials = SignedJwtAssertionCredentials(json_key['client_email'],json_key['private_key'],scope)
gc = gspread.authorize(credentials)

### GOOD DAYS for SEASON 2 ###
if opts.pagenum == 1 and opts.season == 2:
    wkst = gc.open(opts.wkst).sheet1
    if opts.pol == 'xx':
        good_days_epoch1 = wkst.range('C2:C55')
        good_days_epoch2 = wkst.range('C56:C116')
        good_days_epoch3 = wkst.range('C117:C246')
    if opts.pol == 'yy':
        good_days_epoch1 = wkst.range('D2:D55')
        good_days_epoch2 = wkst.range('D56:D116')
        good_days_epoch3 = wkst.range('D117:D246')

    file1 = open('good_days_epoch1_' + opts.pol + '.txt','w')
    file2 = open('good_days_epoch2_' + opts.pol + '.txt','w')
    file3 = open('good_days_epoch3_' + opts.pol + '.txt','w')

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

### GOOD DAYS for SEASON 1 ###
if opts.pagenum == 1 and opts.season == 1:
    wkst = gc.open(opts.wkst).sheet1
    if opts.pol == 'xx':
        good_days_epoch1 = wkst.range('N2:N62')
        good_days_epoch2 = wkst.range('N63:N115')
   
    file1 = open('good_days_epoch1_'+opts.pol+'.txt','w')
    file2 = open('good_days_epoch2_'+opts.pol+'.txt','w')
   
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

### BAD ANTS for SEASON 2 ###
if opts.pagenum == 2 and opts.season == 2:
    wkst = gc.open(opts.wkst).worksheet("Antennas")
    ants = wkst.range('A2:A133')
    epoch2 = wkst.range('H2:H133')
    epoch3 = wkst.range('M2:M133')
    all = wkst.range('B2:B133')

    bad_2 = []
    bad_3 = []
    good_2 = []
    good_3 = []
    bad_all = []
    good_all = []
    for a in range(len(epoch2)):
        val2 = epoch2[a].value
        val3 = epoch3[a].value
        valall = all[a].value
        if val2 == '1':
            bad_2.append(ants[a].value)
        elif val2 == '0':
            good_2.append(ants[a].value)
        if val3 == '1':
            bad_3.append(ants[a].value)
        elif val3 == '0':
            good_3.append(ants[a].value)
        if valall == '1':
            bad_all.append(ants[a].value)
        if valall == '0':
            good_all.append(ants[a].value)
    print 'Epoch 2 Good Antennas:'
    print '   ['+','.join(good_2)+']'
    print 'Epoch 2 Bad Antennas:'
    print '   ['+','.join(bad_2)+']' 
    print 'Epoch 3 Good Antennas:'
    print '   ['+','.join(good_3)+']'
    print 'Epoch 3 Bad Antennas:'
    print '   ['+','.join(bad_3)+']'  
    print 'Cumulative Good Antennas:'
    print '   ['+','.join(good_all)+']'
    print 'Cumulative Bad Antennas:'
    print '   ['+','.join(bad_all)+']'
