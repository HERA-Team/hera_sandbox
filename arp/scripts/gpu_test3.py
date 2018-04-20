# attempt to do matrix multiplication for complex numbers
import pycuda.autoinit
from pycuda import compiler, gpuarray
import numpy as np
import math, time

gpu_template = """
#include <cuComplex.h>

__global__ void MatSqrCmplx(cuFloatComplex *A, cuFloatComplex *C)
{
      const uint A_height = %(A_SHAPE0)s;
      const uint A_width  = %(A_SHAPE1)s;
      const uint row = blockIdx.x * blockDim.x + threadIdx.x;
      const uint col = blockIdx.y * blockDim.y + threadIdx.y;
      cuFloatComplex Csum = make_cuFloatComplex(0,0);

      if (row >= A_height || col >= A_height) return;

      for (int i=0; i < A_width; i++) {
          Csum = cuCaddf(Csum, cuCmulf(A[row * A_width + i], cuConjf(A[col * A_width +i])));
      }
      C[row * A_height + col] = Csum;
}

__global__ void MatMulCmplx(cuFloatComplex *A, cuFloatComplex *B, cuFloatComplex *C)
{
      const uint A_height = %(A_SHAPE0)s;
      const uint A_width  = %(A_SHAPE1)s;
      const uint B_width  = %(B_SHAPE1)s;
      const uint row = blockIdx.x * blockDim.x + threadIdx.x;
      const uint col = blockIdx.y * blockDim.y + threadIdx.y;
      cuFloatComplex Csum = make_cuFloatComplex(0,0);

      if (row >= A_height || col >= B_width) return;

      for (int i=0; i < A_width; i++) {
          Csum = cuCaddf(Csum, cuCmulf(A[row * A_width + i], B[i * B_width + col]));
      }
      C[row * B_width + col] = Csum;
}

__global__ void MatMulFloat(float *A, float *B, float *C)
{
      const uint A_height = %(A_SHAPE0)s;
      const uint A_width  = %(A_SHAPE1)s;
      const uint B_width  = %(B_SHAPE1)s;
      const uint row = blockIdx.x * blockDim.x + threadIdx.x;
      const uint col = blockIdx.y * blockDim.y + threadIdx.y;
      float Csum = 0;

      if (row >= A_height || col >= B_width) return;

      for (int i=0; i < A_width; i++) {
          Csum += A[row * A_width + i] * B[i * B_width + col];
      }
      C[row * B_width + col] = Csum;
}
"""

a_cpu = np.zeros(shape=(32,3145728), dtype=np.complex64)
b_cpu = np.zeros(shape=(3145728,64), dtype=np.complex64)

a_cpu[:,:] = 1+3j
b_cpu[:,:] = 1 + 2j

t1 = time.time()
aa_cpu = np.dot(a_cpu, a_cpu.T.conj())
ab_cpu = np.dot(a_cpu, b_cpu)
cd_cpu = np.dot(a_cpu.real, b_cpu.real)
t_cpu = time.time()-t1

block = (32, 32, 1) # make (32,16,1) for smaller GPUs
a_gpu = gpuarray.to_gpu(a_cpu)
b_gpu = gpuarray.to_gpu(b_cpu)
aa_gpu = gpuarray.empty((a_cpu.shape[0], a_cpu.shape[0]), np.complex64)
ab_gpu = gpuarray.empty((a_cpu.shape[0], b_cpu.shape[1]), np.complex64)
c_gpu = gpuarray.to_gpu(a_cpu.real)
d_gpu = gpuarray.to_gpu(b_cpu.real)
cd_gpu = gpuarray.empty((a_cpu.shape[0], b_cpu.shape[1]), np.float32)

#kernel_code = matmulcmplx_ker % {
#kernel_code = matsqrcmplx_ker % {
gpu_code = gpu_template % {
        'A_SHAPE0': a_cpu.shape[0],
        'A_SHAPE1': a_cpu.shape[1],
        'B_SHAPE1': b_cpu.shape[1],
        }

gpu_module = compiler.SourceModule(gpu_code)
matsqrcmplx = gpu_module.get_function("MatSqrCmplx")
matmulcmplx = gpu_module.get_function("MatMulCmplx")
matmulfloat = gpu_module.get_function("MatMulFloat")
grid = (int(math.ceil(a_cpu.shape[0]/float(block[0]))),int(math.ceil(a_cpu.shape[0]/float(block[1]))))
print grid, block
t1 = time.time()
matsqrcmplx(a_gpu, aa_gpu, grid=grid, block=block)
grid = (int(math.ceil(a_cpu.shape[0]/float(block[0]))),int(math.ceil(b_cpu.shape[1]/float(block[1]))))
print grid, block
matmulcmplx(a_gpu, b_gpu, ab_gpu, grid=grid, block=block)
matmulfloat(c_gpu, d_gpu, cd_gpu, grid=grid, block=block)
t_gpu = time.time()-t1

print("CPU Time ", t_cpu)
print("GPU Time ", t_gpu)

print np.allclose(aa_cpu, aa_gpu.get() ) 
print np.allclose(ab_cpu, ab_gpu.get() ) 
print np.allclose(cd_cpu, cd_gpu.get() ) 
