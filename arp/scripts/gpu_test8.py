import pycuda.autoinit
from pycuda import compiler, gpuarray, driver, cumath
import skcuda.cublas
import numpy as np
import math, time
import aipy

GPU = True
NANT = 128
#NANT = 32
#NANT = 1024
START_JD = 2458000
END_JD = 2458001
INT_TIME = 2160
NSIDE = 512
#NSIDE = 8
NPIX = 12*NSIDE**2

gpu_template = """
#include <cuComplex.h>

__global__ void OneBeamMeasEq(float *A, float *I, float *tau, cuFloatComplex *v)
{
    const uint nant = %(NANT)s;
    const uint npix = %(NPIX)s;
    //const uint tx = threadIdx.x;
    //const uint ty = threadIdx.y;
    const uint ty = threadIdx.x; // switched to make first dim px
    const uint tx = threadIdx.y; // switched to make first dim px
    //const uint row = blockIdx.x * blockDim.x + threadIdx.x;
    //const uint col = blockIdx.y * blockDim.y + threadIdx.y;
    const uint row = blockIdx.y * blockDim.y + threadIdx.y; // switched to make first dim px
    const uint col = blockIdx.x * blockDim.x + threadIdx.x; // switched to make first dim px
    const float freq = %(FREQ)s;
    float amp, phs;
    const uint by = %(BLOCK_Y)s;
    __shared__ float ai[by];

    if (row >= nant || col >= npix) return;
    if (tx == 0)
        ai[ty] = A[col] * I[col];
    __syncthreads(); // make sure all memory is loaded before computing

    amp = ai[ty];
    phs = tau[row,col] * freq;
    v[row,col] = make_cuFloatComplex(amp * cos(phs), amp * sin(phs));
    __syncthreads(); // make everyone used mem before kicking out
}
"""

#block = (1024/NANT, NANT, 1) # make (32,16,1) for smaller GPUs
block = (1024/NANT, NANT, 1) # make (32,16,1) for smaller GPUs
gpu_code = gpu_template % {
        'NANT': NANT,
        'NPIX': NPIX,
        'FREQ': .150,
        'BLOCK_Y': block[1],
        }
gpu_module = compiler.SourceModule(gpu_code)
meq = gpu_module.get_function("OneBeamMeasEq")
#grid = (int(math.ceil(NANT/float(block[0]))),int(math.ceil(NPIX/float(block[1]))))
grid = (int(math.ceil(NPIX/float(block[1]))), int(math.ceil(NANT/float(block[0]))))

times = np.arange(START_JD, END_JD, INT_TIME / aipy.const.s_per_day)

antpos = np.zeros(shape=(NANT,3), dtype=np.float32) # multiply -2pi/c into here
crdtop = np.zeros(shape=(3,NPIX), dtype=np.float32)
A = np.zeros(shape=(1,NPIX), dtype=np.float32)
Isqrt = np.zeros(shape=(1,NPIX), dtype=np.float32)
vis = np.empty(shape=(NANT,NANT), dtype=np.complex64)

antpos_gpu = gpuarray.to_gpu(antpos)
Isqrt_gpu = gpuarray.to_gpu(Isqrt)
crdtop_gpu = gpuarray.empty(shape=crdtop.shape, dtype=crdtop.dtype)
A_gpu = gpuarray.empty(shape=A.shape, dtype=crdtop.dtype)
tau_gpu = gpuarray.empty(shape=(NANT,NPIX), dtype=np.float32)
v_gpu = gpuarray.empty(shape=(NANT,NPIX), dtype=np.complex64)
vis_gpu = gpuarray.empty(shape=vis.shape, dtype=vis.dtype)
h = skcuda.cublas.cublasCreate()
#tau_gpu.fill(0)
#v_gpu.fill(0)
#vis_gpu.fill(0)

print '# Antennas:', NANT
print 'NSIDE:', NSIDE
print 'Starting', time.time()
print grid, block
for ti,jd in enumerate(times):
    print ti,'/',len(times)
    t1 = time.time()
    if GPU:
        A_gpu.set(A)
        crdtop_gpu.set(crdtop)
        t2 = time.time()
        skcuda.cublas.cublasSgemm(h, 'n', 'n', antpos_gpu.shape[0], crdtop_gpu.shape[1], 3, 1., antpos_gpu.gpudata, antpos_gpu.shape[0], crdtop_gpu.gpudata, crdtop_gpu.shape[0], 0., tau_gpu.gpudata, tau_gpu.shape[0])
        meq(A_gpu, Isqrt_gpu, tau_gpu, v_gpu, grid=grid, block=block)
        skcuda.cublas.cublasCgemm(h, 'n', 't', v_gpu.shape[0], v_gpu.shape[0], v_gpu.shape[1], 1., v_gpu.gpudata, v_gpu.shape[0], v_gpu.conj().gpudata, v_gpu.shape[0], 0., vis_gpu.gpudata, vis_gpu.shape[0])
        t3 = time.time()
        vis_gpu.get(vis)
        print t2-t1,t3-t2
        #print aa_cpu
        #print np.allclose(aa_cpu, np.dot(a_cpu, a_cpu.T.conj()))
    else:
        aa_cpu = np.dot(a_cpu, a_cpu.T.conj())
    print time.time() - t1

skcuda.cublas.cublasDestroy(h)
print 'Done', time.time()
