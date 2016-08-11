import unittest
from capo.eor_results import (
        get_pk_from_npz as load_pk, get_k3pk_from_npz as load_k3pk,
        consolidate_bootstraps
        )
from capo.cosmo_units import f212z, c
from capo import pspec
import numpy as np
import os

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
                np.allclose([ref_pk,ref_pk_err], [out_pk,out_pk_err],),
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
        ref_k_mag = np.sqrt( ref_k_par**2 + ref_k_perp**2)
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
                np.allclose([ref_k3pk,ref_k3pk_err], [out_k3pk,out_k3pk_err],),
                msg='Expected pk, pk_err {0} +/- {1} but ' \
                        'retruned pk,pk_err {2} +/- {3}'\
                        ''.format(ref_k3pk,ref_k3pk_err,out_k3pk,out_k3pk_err)
                        )
class Test_bootstrapper(unittest.TestCase):

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))


    def tearDown(self):
        self.path =''

    def test_none_input(self):
        with self.assertRaises(TypeError):
            consolidate_bootstraps()

    def test_blank_list(self):
        with self.assertRaises(TypeError):
            consolidate_bootstraps([])

    def test_empty_string_input(self):
        with self.assertRaises(TypeError):
            consolidate_bootstraps('')

    def test_num_boots(self):
        test_files= [os.path.join(self.path,'test_data/inject_test1')]
        # /test_boot{0:02d}'.format(n)) for n in range(5)]
        ref_boot = 31
        out_dict = consolidate_bootstraps( test_files,
                        save=False,NBOOT=ref_boot,inject=True)
        out_boot = np.shape(out_dict['pCs'])[0]
        self.assertEqual(ref_boot, out_boot)

    def test_load_final_boot(self):
        ref_freq=0.11970443349753694
        test_files= [os.path.join(self.path,'test_data/inject_test1')]
        out_dict = consolidate_bootstraps( test_files,
                        save=False,NBOOT=31,inject=True)
        out_freq = out_dict['freq']
        self.assertEqual(out_freq,ref_freq)

    def test_num_ks(self):
        test_files= [os.path.join(self.path,'test_data/inject_test1')]
        ref_ks = 21
        ref_boot=15
        out_dict = consolidate_bootstraps( test_files,
                    save=False,NBOOT=ref_boot,inject=True)
        out_ks = np.shape(out_dict['pCs'])[1]
        self.assertEqual(out_ks, ref_ks)

    def test_num_injs(self):
        test_files= [os.path.join(self.path,
                'test_data/inject_test{0:01d}'.format(n+1)) for n in range(2) ]
        ref_injs = 2
        ref_boot=15
        out_dict = consolidate_bootstraps( test_files,inject=True,
                    save=False,NBOOT=ref_boot)
        out_injs = np.shape(out_dict['pCs'])[2]
        self.assertEqual(ref_injs, out_injs)

if __name__== '__main__':
    unittest.main()
