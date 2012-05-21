#! /usr/bin/env python
"""
Pulls out the grid 
"""
import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('--correct', action='store_true',
    help='Apply the grid rewire correction on the fly (for use on .uv files)')
opts,args = o.parse_args(sys.argv[1:])
grid_num = {'A0':49,'B0':10,'C0': 9,'D0':22,
            'A1':41,'B1':03,'C1':58,'D1':61,
            'A2':47,'B2':25,'C2':01,'D2':35,
            'A3':19,'B3':48,'C3':04,'D3':18,
            'A4':29,'B4':24,'C4':17,'D4':05,
            'A5':28,'B5':55,'C5':13,'D5':32,
            'A6':34,'B6':27,'C6':56,'D6':30,
            'A7':51,'B7':57,'C7':59,'D7':23  
}
rewire = { 0:49,1:10,2:9,3:22,4:29,5:24,6:17,7:5,
           8:47,9:25,10:1,11:35,12:34,13:27,14:56,15:30,
           16:41,17:3,18:58,19:61,20:28,21:55,22:13,23:32,
           24:19,25:48,26:4,27:18,28:51,29:57,30:59,31:23}
inv_rewire = dict(zip(rewire.values(),rewire.keys()))         
rows = ['A','B','C','D']
cols = range(8)
if opts.correct:
    for k in grid_num:
        grid_num[k] = inv_rewire[grid_num[k]]
#get the EW 16lambda spacing. A0_A1,B0_B1,.. etc
#and the autos
ants = ''
EWs = 0
Ds = 0
As = 0
Dsr  = 0
for col in cols[:-1]:
    for r,row in enumerate(rows):
        ants += "%d_%d,"%(grid_num[str(row)+str(col)],grid_num[str(row)+str(col+1)])
        EWs +=1
        if r<(len(rows)-1):
            print row+str(col)+'_'+rows[r+1]+str(col+1)
            ants += "%d_%d,"%(grid_num[row+str(col)],grid_num[rows[r+1]+str(col+1)])
            Ds += 1
        if r>0:
            ants += "%d_%d,"%(grid_num[row+str(col)],grid_num[rows[r-1]+str(col+1)])
            Dsr +=1

for col in cols:
    for row in rows:
        ants += "%d_%d,"%(grid_num[str(row)+str(col)],grid_num[str(row)+str(col)])
        As += 1
print "ants = ",
print "found %d EW, %d diags (%d down, %d up), %d autos"%(EWs,Ds+Dsr,Ds,Dsr,As)
print ants
for filename in args:
    outfile = os.path.basename(filename)+'G'
    print filename, '->', outfile
    uvi = a.miriad.UV(filename)
    if os.path.exists(outfile):
        print '    File exists... skipping.'
        continue
    a.scripting.uv_selector(uvi, ants=ants)
    #if not opts.pol is None: a.scripting.uv_selector(uvi, pol_str=opts.pol)
    #get all the polarizations!
    for pol in ['xx','yy']:
        a.scripting.uv_selector(uvi, pol_str=pol)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,append2hist='PULL gridspacing:'+' '.join(sys.argv)+'\n')
