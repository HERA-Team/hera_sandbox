import pycuda.autoinit
from pycuda import compiler, gpuarray, driver
import skcuda.cublas
import numpy as np
import math, time
import aipy

GPU = True
NANT = 128
#NANT = 1024
START_JD = 2458000
END_JD = 2458001
INT_TIME = 21600
NSIDE = 512
#NSIDE = 8

times = np.arange(START_JD, END_JD, INT_TIME / aipy.const.s_per_day)

a_cpu = np.zeros(shape=(NANT,12*NSIDE**2), dtype=np.complex64); a_cpu[:,:] = 1.+1j
aa_cpu = np.empty(shape=(NANT,NANT), dtype=a_cpu.dtype)

a_gpu = gpuarray.empty(a_cpu.shape, a_cpu.dtype)
aa_gpu = gpuarray.empty((a_cpu.shape[0], a_cpu.shape[0]), a_cpu.dtype)
h = skcuda.cublas.cublasCreate()

print '# Antennas:', NANT
print 'NSIDE:', NSIDE
print 'Starting', time.time()
for ti,jd in enumerate(times):
    print ti,'/',len(times)
    t1 = time.time()
    if GPU:
        a_gpu.set(a_cpu)
        t2 = time.time()
        #skcuda.cublas.cublasSgemm(h, 'n', 't', a_gpu.shape[0], a_gpu.shape[0], a_gpu.shape[1], 1., a_gpu.gpudata, a_gpu.shape[0], a_gpu.gpudata, a_gpu.shape[0], 0., aa_gpu.gpudata, aa_gpu.shape[0])
        skcuda.cublas.cublasCgemm(h, 'n', 't', a_gpu.shape[0], a_gpu.shape[0], a_gpu.shape[1], 1., a_gpu.gpudata, a_gpu.shape[0], a_gpu.conj().gpudata, a_gpu.shape[0], 0., aa_gpu.gpudata, aa_gpu.shape[0])
        t3 = time.time()
        #driver.memcpy_dtoh(aa_cpu, aa_gpu)
        aa_cpu = aa_gpu.get(aa_cpu)
        print t2-t1,t3-t2,t4-t3
        #print aa_cpu
        #print np.allclose(aa_cpu, np.dot(a_cpu, a_cpu.T.conj()))
    else:
        aa_cpu = np.dot(a_cpu, a_cpu.T.conj())
    print time.time() - t1

skcuda.cublas.cublasDestroy(h)
print 'Done', time.time()
