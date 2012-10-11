import numpy as np
import aipy as a

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
