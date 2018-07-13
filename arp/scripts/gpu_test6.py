import pycuda.autoinit
from pycuda import compiler, gpuarray, driver
import numpy as np
import math, time
import aipy

gpu_template = """
__global__ void MatSqrFloat(float *A, float*C)
{
    const uint A_height = %(A_SHAPE0)s;
    const uint A_width  = %(A_SHAPE1)s;
    const uint tx = threadIdx.x;
    const uint ty = threadIdx.y;
    const uint bx = %(BLOCK_X)s;
    const uint by = %(BLOCK_Y)s;
    const uint smem_depth = %(SMEM_DEPTH)s;
    const uint row = blockIdx.x * blockDim.x + threadIdx.x;
    const uint col = blockIdx.y * blockDim.y + threadIdx.y;
    const uint nchunks = A_width / smem_depth + 1;
    uint chunk_size;
    float Csum = 0;
    __shared__ float a[bx][smem_depth], b[smem_depth][by];

    for (int j=0; j < nchunks; j++) {
        chunk_size = min(smem_depth, A_width - j*smem_depth);
        if (row < A_height) {
            for (int i=ty; i < chunk_size; i += by)
                a[tx][i] = A[row * A_width + j*smem_depth + i];
        }
        if (col < A_height) {
            for (int i=tx; i < chunk_size; i += bx)
                b[i][ty] = A[col * A_width + j*smem_depth + i];
        }
        __syncthreads(); // make sure all memory is loaded before computing

        if (row < A_height && col < A_height) {
            for (int i=0; i < chunk_size; i++)
                Csum += a[tx][i] * b[i][ty];
        }
        __syncthreads(); // make sure everyone's done before loading next chunk
    }
    if (row < A_height && col < A_height)
        C[row * A_height + col] = Csum;
}
"""

GPU = True
NANT = 32
START_JD = 2458000
END_JD = 2458001
INT_TIME = 21600
NSIDE = 512
#NSIDE = 8

times = np.arange(START_JD, END_JD, INT_TIME / aipy.const.s_per_day)

a_cpu = np.zeros(shape=(NANT,12*NSIDE**2), dtype=np.float32); a_cpu[:,:] = 1.
aa_cpu = np.empty(shape=(NANT,NANT), dtype=np.float32)

#block = (32, 32, 1) # make (32,16,1) for smaller GPUs
#block = (2, 512, 1) # make (32,16,1) for smaller GPUs
block = (16, 16, 1) # make (32,16,1) for smaller GPUs
#block = (8, 8, 1) # make (32,16,1) for smaller GPUs
#a_gpu = driver.mem_alloc(a_cpu.nbytes)
#aa_gpu = driver.mem_alloc(aa_cpu.nbytes)
a_gpu = gpuarray.empty(a_cpu.shape, a_cpu.dtype)
aa_gpu = gpuarray.empty((a_cpu.shape[0], a_cpu.shape[0]), a_cpu.dtype)
gpu_code = gpu_template % {
        'A_SHAPE0': a_cpu.shape[0],
        'A_SHAPE1': a_cpu.shape[1],
        'BLOCK_X': block[0],
        'BLOCK_Y': block[1],
        'SMEM_DEPTH': 2**16 / (block[0] + block[1]) / 4 / 2,
        }
gpu_module = compiler.SourceModule(gpu_code)
matsqrcmplx = gpu_module.get_function("MatSqrFloat")
grid = (int(math.ceil(a_cpu.shape[0]/float(block[0]))),int(math.ceil(a_cpu.shape[0]/float(block[1]))))
print grid, block

import time
print '# Antennas:', NANT
print 'NSIDE:', NSIDE
print 'Starting', time.time()
for ti,jd in enumerate(times):
    print ti,'/',len(times)
    if GPU:
        t1 = time.time()
        #driver.memcpy_htod(a_gpu,a_cpu)
        a_gpu.set(a_cpu)
        t2 = time.time()
        matsqrcmplx(a_gpu, aa_gpu, grid=grid, block=block)
        t3 = time.time()
        #driver.memcpy_dtoh(aa_cpu, aa_gpu)
        aa_cpu = aa_gpu.get(aa_cpu)
        t4 = time.time()
        print t2-t1,t3-t2,t4-t3
        #print aa_cpu
        #print np.allclose(aa_cpu, np.dot(a_cpu, a_cpu.T.conj()))
    else:
        aa_cpu = np.dot(a_cpu, a_cpu.T.conj())

print 'Done', time.time()
