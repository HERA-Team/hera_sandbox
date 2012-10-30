import numpy as np
import aipy as a
import ephem
from scipy.special import fdtri as Finv

def bit_flip(f): return np.where(f==1,0,1)

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij
def grid2ij(GRID):
    bls,conj = {},{}
    for row_i in range(GRID.shape[1]):
        for col_i in range(GRID.shape[0]):
            for row_j in range(GRID.shape[1]):
                for col_j in range(GRID.shape[0]):
                    if row_i >= row_j and col_i == col_j: continue
                    sep = '%d,%d' % (col_j-col_i, row_j-row_i)
                    bls[sep] = bls.get(sep,[]) + [(GRID[col_i,row_i],GRID[col_j,row_j])]
    for sep in bls.keys():
        if sep == '0,0' or len(bls[sep]) < 2: del(bls[sep])
    for sep in bls:
        conj[sep] = [i>j for i,j in bls[sep]]
    
    bl_str,bl_conj = {},{}
    for sep in bls:
        bl_str[sep],bl_list = [],[]
        for (i,j),c in zip(bls[sep],conj[sep]):
            if c: i,j = j,i
            bl_list.append(ij2bl(i,j))
            bl_str[sep].append('%d_%d'%(i,j))
            bl_conj[ij2bl(i,j)] = c
        bls[sep] = bl_list
        bl_str[sep] = ','.join(bl_str[sep])
    return bl_str,bl_conj

rad_per_s = 2.*np.pi / a.const.sidereal_day

def lst2bin(lst,bin_width=30.):  
    bin_width *= rad_per_s 
    return np.round(lst / bin_width)
def bin2lst(bin,bin_width=30.):
    return bin / (bin_width * rad_per_s)

def which_lst_files(lst_in,FileDict,tfile=600.): 
    lstout = lst_in + tfile * (2.*np.pi / a.const.sidereal_day)
    outfiles = ['sum','sumsq','wgt']
    _ts = {'sum':[],'sumsq':[],'wgt':[]}
    for i,file in enumerate(FileDict['sum']):
        RAin,RAout = map(ephem.hours,(file.split('.')[-2]).split('_'))
        if lstout < RAin or lst_in > RAout: continue
        if not file in _ts['sum']: 
            for o in outfiles: _ts[o].append(TS[o][i]) 
    return _ts

def flag_F(_d,d1,d2,w,f,alpha = 0.99):
    mean = d1/w
    Svar = (d2/w) - np.abs(mean)**2
    _F = (w/(w+1.))*(np.abs(_d - mean)**2/Svar)
    _Flim = np.array([Finv(2,int(2*(_w-1)),alpha) for _w in w])     
    return _F <= _Flim
