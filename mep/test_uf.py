import numpy as n
import useful_functions as uf
import unittest

class TestPinv(unittest.TestCase):

    def test_projection_matrix(self):
        M = n.diag(n.arange(1,9))
        Z = n.zeros((8,2))
        Z[0,0] = 1; Z[1,1]=1
        PP = uf.projection_matrix(Z)
        #print M.shape, PP.shape
        M_tilde = n.dot(PP.H,n.dot(M,PP))
        self.assertEqual(M_tilde[0,0],0)
        self.assertEqual(M_tilde[1,1],0)

    def test_invertible(self):
        M = n.diag(n.arange(1,9,dtype='float'))
        Minv = uf.pseudo_inverse(M,num_remov=0)
        self.assertTrue(n.all(Minv==n.where(M>0,1/M,0)))

        Minv = uf.pseudo_inverse(M,num_remov=2)
        #print 'remov 2 ', Minv
        self.assertTrue(n.all(Minv==n.diag(n.array([0,0,1/3.,1/4.,1/5.,1/6.,1/7.,1/8.]))))

    def test_pseudo_inverse(self):
        M = n.diag(n.arange(1,9))
        Z = n.zeros((8,2))
        Z[0,0] = 1; Z[1,1]=1
        PP = uf.projection_matrix(Z)
        M_tilde = n.dot(PP.H,n.dot(M,PP))
        Minv = uf.pseudo_inverse(M,num_remov=2)
        M_id = n.dot(Minv,M_tilde)
        proj_id = n.dot(PP.H,n.dot(n.identity(8),PP))
        self.assertTrue(n.all(M_id==proj_id))

    def test_pseudo_remov(self):
        #M = n.diag(n.array([100,10,1,0.1,0.01,0.001]))
        M = n.diag(n.array([0.00001,0.0001,0.001,0.01,0.1,1])[::-1])
        Minv = uf.pseudo_inverse(M,num_remov=None)
        print M
        print Minv 
        Minvp = uf.pseudo_inverse(M,num_remov=2)
        print Minvp
        self.assertTrue(n.allclose(Minv,Minvp))



if __name__=='__main__':
    unittest.main()