import aipy as a, numpy as n, pylab as p
import capo as C
import matplotlib as mpl
import healpy as hp
import useful_functions as uf
import global_sky_model as gsm
import sph_harm_coeffs as shc 
import Q_gsm_error_analysis as qgea 

"""
The monte carlo mpi code writes out a bunch of matrices. Each matrix has rows that 
are vectors of visibilities, so the number of columns is the number of baselines. 
Each row is a different set of visibilities, and each file has 10,000 rows.
"""

def 



