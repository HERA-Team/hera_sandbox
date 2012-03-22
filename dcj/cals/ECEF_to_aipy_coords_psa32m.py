#ipython script to check the positions of the grid
import numpy as n,aipy as a

def length(A):
    return n.sqrt(n.dot(A,A))

F = open('psa32m-20111130.txt')
A = F.readlines()
east_edges = []
west_edges = []
phase_centers  = []
cols = range(8)
rows = ['a','b','c','d']
for col in cols:
    phase_centers.append('a'+str(col)+'c')
    for row in rows:
        east_edges.append(row+str(col) + 'e')
        west_edges.append(row+str(col)+'w')

cpos = [] 
wpos = {}
epos = {}
Names = []
for L in A[1:]:
    name = L.split(',')[0][1:-1]
    try:
        x = float(L.split(',')[3])
        y = float(L.split(',')[4])
        z = float(L.split(',')[5])
    except(ValueError):continue
    #print name
    if name in phase_centers:
        Names.append(name)
        cpos.append([x,y,z])
    if name in east_edges:
        epos[name] = n.array([x,y,z])
    if name in west_edges:
        wpos[name] = n.array([x,y,z])
namesort = n.argsort(Names)

pos = {}
for row in rows:
    for col in cols:
        pos[row+str(col)] = (epos[row+str(col)+'e'] + wpos[row+str(col)+'w'])/2


#long and lat for PSA (radians)
lon = -0.373994485068
lat = -0.536191810965

m= a.coord.rot_m(-1.*lon,n.array([0,0,1]))
m_per_ns = 0.299792458
#antposr = n.dot(m,antpos.transpose()).transpose()

i = 0
print '{'
for row in rows:
    for col in cols:
#        print row+str(col),pos[row+str(col)] - pos['a0']
        #V = (pos[row+str(col)] - pos['a0'])/m_per_ns
        V = n.dot(m,(pos[row+str(col)] - pos['a0']))/m_per_ns
        print i,':',list(V),','
        i +=1
     #   print '----'
print '}'
cpos = n.array(cpos)[namesort].astype(float)
colseps = {}
for col in cols[:-1]:
    colseps[col] = []
    for row in rows:
        colseps[col].append(length(pos[row+str(col)] - pos[row+str(col+1)]))
rowseps = {}
for i,row in enumerate(rows[:-1]):
    rowseps[row] = []
    for col in cols:
        rowseps[row].append(length(pos[row+str(col)] - pos[rows[i+1]+str(col)]))
print "row statistics"
for i,row in enumerate(rows[:-1]):
    print row,rows[i+1],n.mean(rowseps[row]),'+/-',n.std(rowseps[row])
print "col statistics"
for col in cols[:-1]:
    print col,col+1,n.mean(colseps[col]),'+/-',n.std(colseps[col])
