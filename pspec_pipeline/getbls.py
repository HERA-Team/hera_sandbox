#!/usr/bin/env python
import aipy as a
import optparse,sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--sep', '-s', action='store',
             help='Separation type.')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, 1024/100e6, .1, 1024)


ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij
def grid2ij(GRID):
    """
    Create a dictionary of baseline strings '%d,%d'%(i,j) and their grid-separations(in battleship notation. Alongside
    is a dictionary of the baselines you'll need to conjugate.

    Input=GRID: np.array of antenna numbers associated with grid position. The index of the antenna number should be the
    grid-position (modulo stupid python indices). If antenna number 3 is in grid position C4, GRID[4,3] = 3.
    """
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

bl_str,bl_conj = grid2ij(aa.ant_layout)
print bl_str[opts.sep]
