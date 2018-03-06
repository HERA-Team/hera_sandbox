# attempt to do matrix multiplication for complex numbers
import pycuda.autoinit
from pycuda import driver, compiler, gpuarray, tools    
import numpy as np
from time import *

kernel_code_template = """
            #include <cuComplex.h>
    __global__ void MatrixMulKernel(cuFloatComplex *A, cuFloatComplex *B, cuFloatComplex *C)
    {
          const uint wA = %(MATRIX_SIZE)s;
          const uint wB = %(MATRIX_SIZE)s;

          // Block index
          const uint bx = blockIdx.x;
          const uint by = blockIdx.y;

          // Thread index
          const uint tx = threadIdx.x;
          const uint ty = threadIdx.y;

          // Index of the first sub-matrix of A processed by the block
          const uint aBegin = wA * %(BLOCK_SIZE)s * by;
          // Index of the last sub-matrix of A processed by the block
          const uint aEnd   = aBegin + wA - 1;
          // Step size used to iterate through the sub-matrices of A
          const uint aStep = %(BLOCK_SIZE)s;

          // Index of the first sub-matrix of B processed by the block
          const int bBegin = %(BLOCK_SIZE)s * bx;
          // Step size used to iterate through the sub-matrcies of B
          const uint bStep = %(BLOCK_SIZE)s * wB;

          // The element of the block sub-matrix that is computed by the thread
          cuFloatComplex Csub = make_cuFloatComplex(0,0);
          // Loop over all the sub-matrices of A and B required to compute the block sub-matrix
          for (int a = aBegin, b = bBegin;
               a <= aEnd;
           a += aStep, b += bStep)
          {
               // Shared memory for the sub-matrix of A
           __shared__ cuFloatComplex As[%(BLOCK_SIZE)s][%(BLOCK_SIZE)s];
           // Shared memory for the sub-matrix of B
           __shared__ cuFloatComplex Bs[%(BLOCK_SIZE)s][%(BLOCK_SIZE)s];

           // Load the matrices from global memory to shared memory;
           // each thread loads one element of each matrix
           As[ty][tx] = make_cuFloatComplex(cuCrealf(A[a + wA*ty + tx]),cuCimagf(A[a + wA*ty + tx]));
           Bs[ty][tx] = make_cuFloatComplex(cuCrealf(B[b + wB*ty + tx]),cuCimagf(B[b + wA*ty + tx]));

           // Synchronize to make sure the matrices are loaded
           __syncthreads();

           // Multiply the two matrcies together
           // each thread computes one element of the block sub-matrix
           for(int k = 0; k < %(BLOCK_SIZE)s; ++k)
           {
                Csub = cuCaddf(Csub,cuCmulf(As[ty][k],Bs[k][tx]));
           } 

           // Synchronize to make sure that the preceding computation
           // is done before loading two new sub-matrices of A and B in the next iteration
           __syncthreads();
         }

         // Write the block sub-matrix to global memory
         // each thread writes one element
         const uint c = wB * %(BLOCK_SIZE)s * by + %(BLOCK_SIZE)s * bx;
         C[c + wB*ty + tx] = make_cuFloatComplex(cuCrealf(Csub), cuCimagf(Csub));
    }
    """

MATRIX_SIZE = 4
TILE_SIZE  = 2
BLOCK_SIZE = TILE_SIZE

a_cpu = np.zeros(shape=(MATRIX_SIZE,MATRIX_SIZE), dtype=np.complex64)
b_cpu = np.zeros(shape=(MATRIX_SIZE,MATRIX_SIZE), dtype=np.complex64)

a_cpu[:,:] = 1
b_cpu[:,:] = 1 + 2j

# compute reference on the CPU to verify GPU computation
t1 = time()
c_cpu = np.dot(a_cpu, b_cpu)
t2 = time()
t_cpu = t2-t1

# transfer host (CPU) memory to device (GPU) memory
a_gpu = gpuarray.to_gpu(a_cpu)
b_gpu = gpuarray.to_gpu(b_cpu)

# create empty gpuarry for the result (C = A * B)
c_gpu = gpuarray.empty((MATRIX_SIZE, MATRIX_SIZE), np.complex64)

# get the kernel code from the template
# by specifying the constant MATRIX_SIZE
kernel_code = kernel_code_template % {
        'MATRIX_SIZE': MATRIX_SIZE,
        'BLOCK_SIZE': BLOCK_SIZE,
        }

# compile the kernel code
mod = compiler.SourceModule(kernel_code)

# get the kernel function from the compiled module
matrixmul = mod.get_function("MatrixMulKernel")

# call the kernel on the card
t1 = time()
matrixmul(a_gpu, b_gpu, c_gpu, grid=(MATRIX_SIZE/TILE_SIZE, MATRIX_SIZE/TILE_SIZE), block=(TILE_SIZE, TILE_SIZE, 1))
t2 = time()
t_gpu = t2-t1

# print the results
print("-" * 80)
print("Matrix A (GPU): ")
print(a_gpu.get())

print("-" * 80)
print("Matrix B (GPU): ")
print(b_gpu.get())

print("-" * 80)
print("Matrix C (GPU): ")
print(c_gpu.get())

print("-" * 80)
print("Matrix C (CPU): ")
print(c_cpu)

print("-" * 80)
print("CPU-GPU Difference: ")
print(c_cpu-c_gpu.get())

print("CPU Time ", t_cpu)
print("GPU Time ", t_gpu)

np.allclose(c_cpu, c_gpu.get() ) 
