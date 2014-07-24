import numpy as n
import useful_functions as uf
import unittest

class TestPinv(unittest.TestCase):

    def test_projection_matrix(self):
        M = n.diag(n.arange(1,9))
        Z = n.zeros((9,2))
        Z[0,0] = 1; Z[1,1]=1
        PP = uf.projection_matrix(Z)
        M_tilde = n.dot(PP.H,n.dot(M,PP))
        self.assertEqual(M_tilde[0,0],0)
        self.assertEqual(M_tilde[1,1],0)

    def test_invertible(self):
        M = n.diag(n.arange(1,9,dtype='float'))
        Minv = uf.pseudo_inverse(M,num_remov=0)
        self.assertTrue(n.all(Minv==n.where(M>0,1/M,0)))

    def test_pseudo_inverse(self):
        M = n.diag(n.arange(1,9))
        Z = n.zeros((9,2))
        Z[0,0] = 1; Z[1,1]=1
        PP = uf.projection_matrix(Z)
        M_tilde = n.dot(PP.H,n.dot(M,PP))
        Minv = uf.pseudo_inverse(M,num_remov=2)
        M_id = n.dot(Minv,M_tilde)
        proj_id = n.dot(PP.H,n.dot(n.identity(8),PP))
        self.assertTrue(n.all(M_id==proj_id))




if __name__=='__main__':
    unittest.main()