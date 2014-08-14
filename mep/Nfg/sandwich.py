#! /usr/bin/env python

import numpy as np, sys

def vect2sq(n,vect):
    sq = np.zeros((n,n),dtype=complex)
    for i in range(n):
        for j in range(n):
            if i >= j:
                sq[i,j] = vect[j+i*(i+1)/2]
            else:
                sq[i,j] = np.conj(vect[i+j*(j+1)/2])
    return sq

if __name__=='__main__':
    nside = int(sys.argv[1])
    npix = 12 * nside * nside
    KvectFname = sys.argv[2]
    GmatrixFname = sys.argv[3]
    outputFname = sys.argv[4]

    Gmatrix = np.load(GmatrixFname)
    Kvect = np.load(KvectFname)
    Kmatrix = vect2sq(npix,Kvect)

    print "Doing first multiplication..."
    GKmatrix = np.einsum('ijm,jk',Gmatrix,Kmatrix)
    print "Doing second multiplication..."
    GKGdagger = np.einsum('ijm,jk',GKmatrix,np.conj(Gmatrix.T))

    np.save(outputFname,GKGdagger)
    for ii,fq in enumerate(fqs):
        n.savez_compressed('{0}_fq_{1}'.format(outputFname,fq),matrix=GKGdagger[:,:,ii])
