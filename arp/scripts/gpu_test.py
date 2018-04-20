import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np

kernels = SourceModule("""
__global__ void custom_kernel(float *g_y, float *g_x)
{
const int i = blockDim.x * blockIdx.x + threadIdx.x;
const float x = g_x[i];
g_y[i] = cos(x)*exp(sin(x)-sqrt(x*x));
}
""")

custom_kernel = kernels.get_function("custom_kernel");
size = 5120000
block_size = 512 # design a 1d block and grid structure
grid_size = size/block_size
block = (block_size,1,1)
grid = (grid_size,1)

X = np.linspace(1,size,size).astype(np.float32)
X_gpu = gpuarray.to_gpu(X)
Y_gpu = gpuarray.empty_like(X_gpu) # 1. transfer to GPU
custom_kernel(Y_gpu, X_gpu, block=block, grid=grid) # 2. execute kernel

answer = gpuarray.sum(Y_gpu) # 3. use PyCUDA sum reduction
print answer.get()

print np.sum(np.cos(X) * np.exp(np.sin(X)-np.sqrt(X*X)))
