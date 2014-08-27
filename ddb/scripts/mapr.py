from subprocess import call

mfp = open('fx_to_ant.dat','r')
mapr = []
for line in mfp:
    d = line.split()
    dd = [d[0],d[2],d[4],d[6],d[8]]
    mapr.append(dd)
mfp.close()

mapf2o = {'A1':0,'A2':1,'A3':2,'A4':3,'B1':4,'B2':5,'B3':6,'B4':7,
          'C1':8,'C2':9,'C3':10,'C4':11,'D1':12,'D2':13,'D3':14,'D4':15,
          'E1':16,'E2':17,'E3':18,'E4':19,'F1':20,'F2':21,'F3':22,'F4':23,
          'G1':24,'G2':25,'G3':26,'G4':27,'H1':28,'H2':29,'H3':30,'H4':31}


def readLevels(name='level.txt'):
    lfp = open(name,'r')
    levels = {}
    for line in lfp:
        d = line.split()
        pf = d[0][1:]
        levels[pf] = d[2:]

def f2r(f,name='refresh'):
    if name=='refresh':
        name = 'levels.tmp'
        getCurrent(name)
    readLevels(name)
    if f[0]=='p':
        del(f[0])
    key = f[0:2]
    n = mapf2o[f[2:4]]
    r = None
    for pf in mapr:
        if pf[0][0:4] == f:
            r = pf[2]
    print f,r,levels[key][n]


def r2f(r,name='refresh'):
    if name=='refresh':
        name = 'levels.tmp'
        getCurrent(name)
    for pf in mapr:
        if pf[2][0:2]==r:
            f = pf[0]
            f2r(f,name)

def getCurrent(name='level.tmp'):
    cmd = 'curl http://paper1:3000/instruments/psa256/levels.txt > '+name
    call(cmd)
