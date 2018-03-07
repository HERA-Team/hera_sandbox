import pycuda.autoinit
from pycuda import compiler, gpuarray, driver, cumath
import skcuda.cublas
import numpy as np
import math, time
import aipy

GPU = True
NANT = 128
#NANT = 8
START_JD = 2458000
END_JD = 2458001
INT_TIME = 21600
NSIDE = 512
#NSIDE = 8
NPIX = 12*NSIDE**2
FREQ = .150
BEAM_PX = 63

gpu_template = """
#include <cuComplex.h>

__device__
inline float lerp(float v0, float v1, float t) {
    return fma(t, v1, fma(-t, v0, v0));
}

// XXX need interpolation
texture<float, cudaTextureType2D, cudaReadModeElementType> Atex;
__shared__ float ai[%(BLOCK_Y)s];

__global__ void InterpolateOneBeam(float *top, float *A)
{
    const uint npix = %(NPIX)s;
    const uint col = blockIdx.x * blockDim.x + threadIdx.x;
    const uint beam_px = %(BEAM_PX)s;
    float bm_x, bm_y, px, py, fx, fy;
    if (col >= npix) return;
    if (top[2*npix+col] > 0) {
        bm_x = beam_px * (0.5 * top[col] + 0.5);
        bm_y = beam_px * (0.5 * top[npix+col] + 0.5);
        bm_x -= 0.5f;
        bm_y -= 0.5f;
        px = floorf(bm_x);   // integer position
        py = floorf(bm_y);
        fx = bm_x - px;    // fractional position
        fy = bm_y - py;
        px += 0.5f;
        py += 0.5f;
        A[col] = lerp(lerp(tex2D(Atex,px,py),      tex2D(Atex,px+1.0f,py),fx),
                 lerp(tex2D(Atex,px,py+1.0f), tex2D(Atex,px+1.0f,py+1.0f),fx), fy);
    } else {
        A[col] = 0;
    }
}

__global__ void OneBeamMeasEq(float *top, float *A, float *I, float *tau, cuFloatComplex *v)
{
    const uint nant = %(NANT)s;
    const uint npix = %(NPIX)s;
    const uint ty = threadIdx.x; // switched to make first dim px
    const uint tx = threadIdx.y; // switched to make second dim ant
    const uint row = blockIdx.y * blockDim.y + threadIdx.y; // switched to make first dim px
    const uint col = blockIdx.x * blockDim.x + threadIdx.x; // switched to make second dim ant
    const uint beam_px = %(BEAM_PX)s;
    const float freq = %(FREQ)s;
    float amp, phs;

    if (row >= nant || col >= npix) return;
    if (tx == 0)
        ai[ty] = A[col] * I[col];

    __syncthreads(); // make sure all memory is loaded before computing

    amp = ai[ty];
    phs = tau[row*npix + col] * freq;
    v[row*npix + col] = make_cuFloatComplex(amp * cos(phs), amp * sin(phs));
    __syncthreads(); // make everyone used mem before kicking out
}
"""

#block = (1024/NANT, NANT, 1) # make (32,16,1) for smaller GPUs
block = (1024/NANT, NANT, 1) # make (32,16,1) for smaller GPUs
gpu_code = gpu_template % {
        'NANT': NANT,
        'NPIX': NPIX,
        'FREQ': FREQ,
        'BEAM_PX': BEAM_PX,
        #'BLOCK_Y': block[1],
        'BLOCK_Y': block[0], # switched to make first dim npix
        }
gpu_module = compiler.SourceModule(gpu_code)
beam_interpolate = gpu_module.get_function("InterpolateOneBeam")
meq = gpu_module.get_function("OneBeamMeasEq")
texref = gpu_module.get_texref("Atex")
grid = (int(math.ceil(NPIX/float(block[0]))),int(math.ceil(NANT/float(block[1]))))

times = np.arange(START_JD, END_JD, INT_TIME / aipy.const.s_per_day)


np.random.seed(0)
antpos = np.zeros(shape=(NANT,3), dtype=np.float32) # multiply -2pi/c into here
antpos[:,0] = 1
eq2top = np.array([[1.,0,0],[0,1,0],[0,0,1]], dtype=np.float32)
crdeq = np.random.uniform(size=(3,NPIX)).astype(np.float32)
Atex = np.ones(shape=(BEAM_PX,BEAM_PX), dtype=np.float32) # X is horizontal (2nd dim), Y is vertical (1st dim)
Atex *= 0
#Atex[31,31] = 2
#Atex[0,31] = 2
Atex[31,0] = 2
Isqrt = np.zeros(shape=(1,NPIX), dtype=np.float32)
vis = np.empty(shape=(NANT,NANT), dtype=np.complex64)

driver.matrix_to_texref(Atex, texref, order='C') # never changes
antpos_gpu = gpuarray.to_gpu(antpos) # never changes
Isqrt_gpu = gpuarray.to_gpu(Isqrt) # never changes
A_gpu = gpuarray.empty(shape=(NPIX,), dtype=np.float32)
crdeq_gpu = gpuarray.to_gpu(crdeq) # never changes
eq2top_gpu = gpuarray.to_gpu(eq2top) # sent from CPU each jd
crdtop_gpu = gpuarray.empty(shape=(3,NPIX), dtype=np.float32)
tau_gpu = gpuarray.empty(shape=(NANT,NPIX), dtype=np.float32)
v_gpu = gpuarray.empty(shape=(NANT,NPIX), dtype=np.complex64)
vis_gpu = gpuarray.empty(shape=(NANT,NANT), dtype=np.complex64)
h = skcuda.cublas.cublasCreate()
#tau_gpu.fill(0)
#vis_gpu.fill(0)
#v = np.zeros(shape=(NANT,NPIX), dtype=np.complex64)
#for ant in xrange(NANT):
#    if ant % 2 == 1: v[ant,ant/2] = 1
#    else: v[ant,ant/2] = 1+1j
crdtop = np.zeros((3,NPIX), dtype=np.float32)
crdtop[2] = 1
crdtop[0,0] = -1

print '# Antennas:', NANT
print 'NSIDE:', NSIDE
print 'Starting', time.time()
print grid, block
for ti,jd in enumerate(times):
    print ti,'/',len(times)
    t1 = time.time()
    if GPU:
        eq2top_gpu.set(eq2top)
        t2 = time.time()
        skcuda.cublas.cublasSgemm(h, 'n', 'n', NPIX, 3, 3, 1., crdeq_gpu.gpudata, NPIX, eq2top_gpu.gpudata, 3, 0., crdtop_gpu.gpudata, NPIX)
        #print np.allclose(crdtop_gpu.get(), crdeq)
        skcuda.cublas.cublasSgemm(h, 'n', 'n', NPIX, NANT, 3, 1., crdtop_gpu.gpudata, NPIX, antpos_gpu.gpudata, 3, 0., tau_gpu.gpudata, NPIX)
        #print np.allclose(tau_gpu.get()[0], crdeq[0])
        # XXX still need to test meq
        #crdtop_gpu.set(crdtop)
        beam_interpolate(crdtop_gpu, A_gpu, grid=(NPIX/1024,1), block=(1024,1,1))
        #print A_gpu.get()
        #Isqrt_gpu.fill(1)
        #tau_gpu.fill(0)
        #meq(crdtop_gpu, Isqrt_gpu, tau_gpu, v_gpu, grid=grid, block=block)
        meq(crdtop_gpu, A_gpu, Isqrt_gpu, tau_gpu, v_gpu, grid=grid, block=block)
        #print v_gpu.get()
        #v_gpu.set(v)
        # transpose below incurs about 20% overhead
        skcuda.cublas.cublasCgemm(h, 't', 'n', NANT, NANT, NPIX, 1., v_gpu.conj().gpudata, NPIX, v_gpu.gpudata, NPIX, 0., vis_gpu.gpudata, NANT)
        t3 = time.time()
        vis_gpu.get(vis)
        print t2-t1,t3-t2
        #print vis
        #print np.allclose(aa_cpu, np.dot(a_cpu, a_cpu.T.conj()))
    else:
        aa_cpu = np.dot(a_cpu, a_cpu.T.conj())
    print time.time() - t1

skcuda.cublas.cublasDestroy(h)
print 'Done', time.time()
