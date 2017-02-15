import sys, os, glob

JDs = sys.argv[1:]

pols = ['xx','xy','yx','yy']

for JD in JDs:
    equality = []
    stor = ''
    for pp in pols:
        pp_uvcRRE = sorted(glob.glob('/data4/paper/2014EoR/pot3/%s/*%s.uvcRRE'%(JD,pp)))
        pp_uvcRREc= sorted(glob.glob('/data4/paper/2014EoR/pot3/%s/*%s.uvcRREc'%(JD,pp)))
        if len(pp_uvcRREc)==len(pp_uvcRRE): equality.append(True)
        else: equality.append(False)
        stor+='%i/%i %s RREc:RRE \n'%(len(pp_uvcRREc),len(pp_uvcRRE),pp)
    if not all(equality): 
        print JD
        #print stor
        
        
