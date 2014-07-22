import unittest
import numpy as n
import useful_functions as uf 
import Q_gsm_error_analysis as qgea

class TestRecoverGS(unittest.TestCase):
    def test_recover_diag(self):
        Q = n.diag(n.arange(1,9))
        N = n.identity(8)
        a = n.arange(2,17,2)
        y = n.dot(Q,a)+uf.rand_from_covar(N)
        a,ahat,err = qgea.test_recover_alms(y,Q,N,a,num_remov=-4)
        print "diag"
        print a
        print ahat
        #self.assertTrue(n.all(a==ahat))

    def test_recover_wind(self):
        num=100
        Q = n.diag(n.arange(1,9))
        N = n.identity(8)
        a = n.resize(n.arange(2,17,2),(8,1))
        a1 = n.resize(n.arange(2,17,2),(8,1))
        ahat_avg = 0
        for ii in range(num):
            n_vec = n.resize(uf.rand_from_covar(N),(8,1))
            y = n.dot(Q,a)+n_vec
            self.assertTrue(n.all(y.shape==(8,1)))
            a,ahat,err = qgea.test_recover_alms(y,Q,N,a,num_remov=-4)        
            #print ahat
            ahat_avg += ahat
            self.assertTrue(n.all(y==n.dot(Q,a1)+n_vec))
        ahat_avg = ahat_avg/num
        W = qgea.window_fn_matrix(Q,N,num_remov=-4)
        print "Window"
        print n.dot(W,a)
        print ahat_avg
        self.assertTrue(n.all(n.abs(n.dot(W,a)-ahat_avg)<0.1))


if __name__=='__main__':
    unittest.main()
