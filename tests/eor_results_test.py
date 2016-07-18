import unittest
from capo.eor_results import get_pk_from_npz as load_pk, get_k3pk_from_npz as load_k3pk
from capo.cosmo_units import f212z, c
from capo import pspec
import numpy as n

test_data_dir='test_data/'
test_data_file= test_data_dir+'test_eor_results.npz'

class Test_eor_loader(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_z_from_pk(self):
        ref_z=f212z(119.7e6)
        out_z,_,_,_ = load_pk(test_data_file,verbose=False)
        self.assertEqual(ref_z, out_z,
                msg='Expected z {0} but ' \
                        'retruned z {1}'.format(ref_z,out_z)
                        )

    def test_k_from_pk(self):
        ref_k=.2
        _,out_k,_,_ = load_pk(test_data_file,verbose=False)
        self.assertEqual(ref_k, out_k,
                msg='Expected k_par {0} but ' \
                        'retruned k_par {1}'.format(ref_k,out_k)
                        )

    def test_pk_from_pk(self):
        ref_pk = 36752
        ref_pk_err = 13987
        _,_,out_pk,out_pk_err = load_pk(test_data_file,verbose=False)
        self.assertTrue(
                n.allclose([ref_pk,ref_pk_err], [out_pk,out_pk_err],),
                msg='Expected pk, pk_err {0} +/- {1} but ' \
                        'retruned pk,pk_err {2} +/- {3}'\
                        ''.format(ref_pk,ref_pk_err,out_pk,out_pk_err)
                        )

    def test_z_from_k3pk(self):
        ref_z=f212z(119.7e6)
        out_z,_,_,_ = load_k3pk(test_data_file,verbose=False)
        self.assertEqual(ref_z, out_z,
                msg='Expected z {0} but ' \
                        'retruned z {1}'.format(ref_z,out_z)
                        )

    def test_k_from_k3pk(self):
        ref_k_par=.2
        ref_z=f212z(119.7e6)
        ref_umag = 30/(3e8/(119.7*1e6))
        ref_k_perp = ref_umag*pspec.dk_du(ref_z)
        ref_k_mag = n.sqrt( ref_k_par**2 + ref_k_perp**2)
        _,out_k_mag,_,_ = load_k3pk(test_data_file,verbose=False)
        self.assertEqual(ref_k_mag, out_k_mag,
                msg='Expected k_mag {0} but ' \
                        'retruned k_mag {1}'.format(ref_k_mag,out_k_mag)
                        )

    def test_k3pk_from_k3pk(self):
        ref_k3pk = 992
        ref_k3pk_err = 1003
        _,_,out_k3pk,out_k3pk_err = load_k3pk(test_data_file,verbose=False)
        self.assertTrue(
                n.allclose([ref_k3pk,ref_k3pk_err], [out_k3pk,out_k3pk_err],),
                msg='Expected pk, pk_err {0} +/- {1} but ' \
                        'retruned pk,pk_err {2} +/- {3}'\
                        ''.format(ref_k3pk,ref_k3pk_err,out_k3pk,out_k3pk_err)
                        )

if __name__== '__main__':
    unittest.main()
