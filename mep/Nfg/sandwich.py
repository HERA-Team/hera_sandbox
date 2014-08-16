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
    print "I live!"
    nside = int(sys.argv[1])
    npix = 12 * nside * nside
    KvectFname = sys.argv[2]
    DataFolder = sys.argv[3]
    GmatrixFname = sys.argv[4]
    outputFname = sys.argv[5]

    Gmatrix = np.load('{0}/{1}.npz'.format(DataFolder,GmatrixFname))
    Gmatrix = Gmatrix['arr_0']
    print 'Gmatrix loaded...'
    Kvect = np.load(KvectFname)
    print 'Kvect loaded...'
    Kmatrix = vect2sq(npix,Kvect)

    print "Doing first multiplication..."
    print Gmatrix.shape
    print Kmatrix.shape
    GKmatrix = np.einsum('ijm,jk',Gmatrix,Kmatrix)
    print "Doing second multiplication..."
    print GKmatrix.shape
    GKGdagger = np.einsum('ijm,jk',GKmatrix,np.conj(Gmatrix.T))

#    np.save(outputFname,GKGdagger)
    for ii,fq in enumerate(fqs):
        n.savez_compressed('{0}/{1}_fq_{2}'.format(DataFolder,outputFname,fq),matrix=GKGdagger[:,:,ii])
