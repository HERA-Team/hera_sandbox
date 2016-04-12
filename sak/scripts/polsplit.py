import os, sys, glob

pols = ['xx','xy','yx','yy']

for FILE in sys.argv[1:]:
    print FILE
    
    spFILE = FILE.split('/')
    F = os.path.basename(FILE)
    _polnames=[]
    
    for p in pols:
        spF = F.split('.') 
        spF.insert(3,p)
        jF = '.'.join(spF)
        spFILE[-1] = jF
        jFILE = '/'.join(spFILE)
        print 'pull_antpols.py -p %s %s'%(p,FILE)
        os.system('pull_antpols.py -p %s %s'%(p,FILE))
        print 'mv %sA %s'%(FILE,jFILE) 
        os.system('mv %sA %s'%(FILE,jFILE))
        print 'rm -rf %sA'%FILE
        os.system('rm -rf %sA'%FILE)
