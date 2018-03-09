import pycuda.autoinit
from pycuda import compiler, gpuarray, driver, cumath
import skcuda.cublas
import numpy as np
import math, time
import aipy

NTHREADS = 1024 # make 512 for smaller GPUs
MAX_MEMORY = 2**29
NANT = 256
START_JD = 2458000
END_JD = 2458001
INT_TIME = 21600
NSIDE = 512
NPIX = 12*NSIDE**2 # this is assumed to always be bigger than NTHREADS
CHUNK = max(2,2**int(math.ceil(np.log2(float(NANT*NPIX) / MAX_MEMORY))))
NPIXC = NPIX / CHUNK
FREQ = .150
BEAM_PX = 63

gpu_template = """
#include <cuComplex.h>

// Linearly interpolate between [v0,v1] for t=[0,1]
// v = v0 * (1-t) + v1 * t = t*v1 + (-t*v0 + v0)
// Runs on GPU only
__device__
inline float lerp(float v0, float v1, float t) {
    return fma(t, v1, fma(-t, v0, v0));
}

// Texture storing beam response on (x=sin th_x, y=sin th_y, nant) grid
// for fast lookup from multiple threads.  Suggest setting first 2 dims of
// bm_tex to an odd number to get pixel centered on zenith
texture<float, cudaTextureType3D, cudaReadModeElementType> bm_tex;

// Shared memory for storing per-antenna results to be reused among all ants,
// avoiding a rush on global memory.
//__shared__ float sh_buf[%(BLOCK_Y)s,4];
__shared__ float sh_buf[%(BLOCK_Y)s*5];

// Interpolate bm_tex[x,y] at top=(x,y,z) coordinates and store answer in A
__global__ void InterpolateBeam(float *top, float *A)
{
    const uint nant = %(NANT)s;
    const uint npix = %(NPIX)s;
    const uint ty = threadIdx.x; // switched to make first dim px
    const uint tx = threadIdx.y; // switched to make second dim ant
    const uint col = blockIdx.x * blockDim.x + threadIdx.x;
    const uint ant = blockIdx.y * blockDim.y + threadIdx.y;
    const uint beam_px = %(BEAM_PX)s;
    float bm_x, bm_y, px, py, pz, fx, fy, top_z;
    if (col >= npix || ant >= nant) return;
    if (tx == 0) // buffer top_z for all threads
        sh_buf[ty+%(BLOCK_Y)s * 4] = top[2*npix+col];
    __syncthreads(); // make sure top_z exists for all threads
    top_z = sh_buf[ty+%(BLOCK_Y)s * 4];
    if (tx == 0 && top_z > 0) { // buffer x interpolation for all threads
        bm_x = beam_px * (0.5 * top[col] + 0.5) - 0.5f;
        px = floorf(bm_x);   // integer position
        sh_buf[ty+%(BLOCK_Y)s * 0] = bm_x - px; // fx, fractional position
        sh_buf[ty+%(BLOCK_Y)s * 2] = px + 0.5f; // px, pixel index
    }
    if (tx == 1 && top_z > 0) { // buffer y interpolation for all threads
        bm_y = beam_px * (0.5 * top[npix+col] + 0.5) - 0.5f;
        py = floorf(bm_y);
        sh_buf[ty+%(BLOCK_Y)s * 1] = bm_y - py; // fy, fractional position
        sh_buf[ty+%(BLOCK_Y)s * 3] = py + 0.5f; // py, pixel index
    }
    __syncthreads(); // make interpolation exists for all threads
    if (top_z > 0) {
        fx = sh_buf[ty+%(BLOCK_Y)s * 0];
        fy = sh_buf[ty+%(BLOCK_Y)s * 1];
        px = sh_buf[ty+%(BLOCK_Y)s * 2];
        py = sh_buf[ty+%(BLOCK_Y)s * 3];
        pz = ant + 0.5f;
        A[ant*npix+col] = lerp(lerp(tex3D(bm_tex,px,py,pz),      tex3D(bm_tex,px+1.0f,py,pz),fx),
               lerp(tex3D(bm_tex,px,py+1.0f,pz), tex3D(bm_tex,px+1.0f,py+1.0f,pz),fx), fy);
    } else {
        A[ant*npix+col] = 0;
    }
    __syncthreads(); // make sure everyone used mem before kicking out
}

// Compute A*I*exp(ij*tau*freq) for all antennas, storing output in v
__global__ void MeasEq(float *A, float *I, float *tau, float freq, cuFloatComplex *v)
{
    const uint nant = %(NANT)s;
    const uint npix = %(NPIX)s;
    const uint ty = threadIdx.x; // switched to make first dim px
    const uint tx = threadIdx.y; // switched to make second dim ant
    const uint row = blockIdx.y * blockDim.y + threadIdx.y; // second thread dim is ant
    const uint col = blockIdx.x * blockDim.x + threadIdx.x; // first thread dim is px
    const uint beam_px = %(BEAM_PX)s;
    float amp, phs;

    if (row >= nant || col >= npix) return;
    if (tx == 0)
        sh_buf[ty] = I[col];
    __syncthreads(); // make sure all memory is loaded before computing
    amp = A[row*npix + col] * sh_buf[ty];
    phs = tau[row*npix + col] * freq;
    v[row*npix + col] = make_cuFloatComplex(amp * cos(phs), amp * sin(phs));
    __syncthreads(); // make sure everyone used mem before kicking out
}
"""

def numpy3d_to_array(np_array):
    '''Copy a 3D numpy array into a 3D pycuda array that can be used to set a texture.
    (For some reason, gpuarrays can't be used to do that directly).  A transpose
    happens implicitly, so the resulting array has dimensions (w,h,d).'''
    import pycuda.autoinit
    d, h, w = np_array.shape
    descr = driver.ArrayDescriptor3D()
    descr.width = w
    descr.height = h
    descr.depth = d
    descr.format = driver.dtype_to_array_format(np_array.dtype)
    descr.num_channels = 1
    descr.flags = 0
    device_array = driver.Array(descr)
    copy = driver.Memcpy3D()
    copy.set_src_host(np_array)
    copy.set_dst_array(device_array)
    copy.width_in_bytes = copy.src_pitch = np_array.strides[1]
    copy.src_height = copy.height = h
    copy.depth = d
    copy()
    return device_array

# blocks of threads are mapped to (pixels,ants,freqs)
block = (max(1,NTHREADS/NANT), min(NTHREADS,NANT), 1)
grid = (int(math.ceil(NPIXC/float(block[0]))),int(math.ceil(NANT/float(block[1]))))
gpu_code = gpu_template % {
        'NANT': NANT,
        'NPIX': NPIXC,
        'BEAM_PX': BEAM_PX,
        'BLOCK_Y': block[0], # switched to make first dim npix
        }
gpu_module = compiler.SourceModule(gpu_code)
beam_interpolate = gpu_module.get_function("InterpolateBeam")
meq = gpu_module.get_function("MeasEq")
texref = gpu_module.get_texref("bm_tex")
h = skcuda.cublas.cublasCreate()

# Initialization of values on CPU side
np.random.seed(0)
antpos = np.zeros(shape=(NANT,3), dtype=np.float32) # multiply -2pi/c into here
antpos[:,0] = 1
eq2top = np.array([[1.,0,0],[0,1,0],[0,0,1]], dtype=np.float32)
#crdeq = np.random.uniform(size=(3,NPIX)).astype(np.float32)
crdeq = np.zeros(shape=(3,NPIX), dtype=np.float32)
# note that bm_tex is transposed relative to the cuda texture buffer
bm_tex = np.ones(shape=(NANT,BEAM_PX,BEAM_PX), dtype=np.float32) # X is 3rd dim, Y is 2nd dim
bm_tex *= 0
bm_tex[:,31,0] = 2
Isqrt = np.zeros(shape=(1,NPIX), dtype=np.float32)
Isqrt[0] = 1
#crdtop = np.zeros((3,NPIX), dtype=np.float32)
crdeq[2] = 1
crdeq[0,0] = -1

# Transfer of values to GPU
texref.set_array(numpy3d_to_array(bm_tex)) # never changes, transpose happens in copy so cuda bm_tex is (BEAM_PX,BEAM_PX,NANT)
antpos_gpu = gpuarray.to_gpu(antpos) # never changes, set to -2*pi*antpos/c
#Isqrt_gpu = gpuarray.to_gpu(Isqrt) # never changes
Isqrt_gpu = gpuarray.empty(shape=(1,NPIXC), dtype=np.float32)
A_gpu = gpuarray.empty(shape=(NANT,NPIXC), dtype=np.float32) # will be set on GPU by beam_interpolate
#crdeq_gpu = gpuarray.to_gpu(crdeq) # never changes
crdeq_gpu = gpuarray.empty(shape=(3,NPIXC), dtype=np.float32)
eq2top_gpu = gpuarray.to_gpu(eq2top) # sent from CPU each jd
crdtop_gpu = gpuarray.empty(shape=(3,NPIXC), dtype=np.float32) # will be set on GPU
tau_gpu = gpuarray.empty(shape=(NANT,NPIXC), dtype=np.float32) # will be set on GPU
#v_gpu = gpuarray.empty(shape=(NANT,NPIXC), dtype=np.complex64) # will be set on GPU
v_gpus = [gpuarray.empty(shape=(NANT,NPIXC), dtype=np.complex64) for i in xrange(CHUNK)]
#vis_gpu = gpuarray.empty(shape=(NANT,NANT), dtype=np.complex64) # will be set on GPU and returned

#tau = np.zeros((NANT,NPIXC), dtype=np.float32)
#tau_gpu.fill(0)
#vis_gpu.fill(0)
#v = np.zeros(shape=(NANT,NPIXC), dtype=np.complex64)
#for ant in xrange(NANT):
#    if ant % 2 == 1: v[ant,ant/2] = 1
#    else: v[ant,ant/2] = 1+1j

print '# Antennas:', NANT
print 'NSIDE:', NSIDE
print 'CHUNK:', CHUNK
print 'Starting', time.time()
print grid, block
times = np.arange(START_JD, END_JD, INT_TIME / aipy.const.s_per_day)
event_order = ('start','upload','eq2top','tau','interpolate','meq','vis','end')
#events = [{e:driver.Event() for e in event_order} for i in xrange(CHUNK)]
#streams = [driver.Stream() for i in xrange(CHUNK)]
#for c in xrange(CHUNK):
#    for e in events[c].values(): e.record()
#crdeq_gpu.set_async(crdeq[:,:NPIXC]) # initial for CHUNK=1 case
#Isqrt_gpu.set_async(Isqrt[:,:NPIXC]) # initial for CHUNK=1 case
viss = [np.empty(shape=(NANT,NANT), dtype=np.complex64) for i in xrange(CHUNK)]
# XXX need a separate notion of chunks of data and streams processing it
vis_gpus = [gpuarray.empty(shape=(NANT,NANT), dtype=np.complex64) for i in xrange(CHUNK)]
for ti,jd in enumerate(times):
    print ti,'/',len(times)
    t0 = time.time()
    eq2top_gpu.set(eq2top)
    events = [{e:driver.Event() for e in event_order} for i in xrange(CHUNK)]
    streams = [driver.Stream() for i in xrange(CHUNK)]
    for c in xrange(CHUNK):
        prev_c = (CHUNK-1+c) % CHUNK
        print c,':',CHUNK
        if c > 0: streams[c].wait_for_event(events[prev_c]['meq'])
        events[c]['start'].record(streams[c])
        crdeq_gpu.set_async(crdeq[:,c*NPIXC:(c+1)*NPIXC], stream=streams[c])
        Isqrt_gpu.set_async(Isqrt[:,c*NPIXC:(c+1)*NPIXC], stream=streams[c])
        events[c]['upload'].record(streams[c])
        skcuda.cublas.cublasSetStream(h, streams[c].handle)
        # compute crdtop = dot(eq2top,crdeq)
        # cublas arrays are in Fortran order, so P=M*N is actually peformed as P.T = N.T * M.T
        skcuda.cublas.cublasSgemm(h, 'n', 'n', NPIXC, 3, 3, 1., crdeq_gpu.gpudata, NPIXC, eq2top_gpu.gpudata, 3, 0., crdtop_gpu.gpudata, NPIXC)
        events[c]['eq2top'].record(streams[c])
        ##print np.allclose(crdtop_gpu.get(), crdeq)
        # compute tau = dot(antpos,crdtop)
        skcuda.cublas.cublasSgemm(h, 'n', 'n', NPIXC, NANT, 3, 1., crdtop_gpu.gpudata, NPIXC, antpos_gpu.gpudata, 3, 0., tau_gpu.gpudata, NPIXC)
        events[c]['tau'].record(streams[c])
        ##print np.allclose(tau_gpu.get()[0], crdeq[0])
        # interpolate bm_tex at specified topocentric coords, store interpolation in A
        # threads are parallelized across pixel axis
        beam_interpolate(crdtop_gpu, A_gpu, grid=grid, block=block, stream=streams[c])
        events[c]['interpolate'].record(streams[c])
        # compute v = A * I * exp(1j*tau*freq)
        meq(A_gpu, Isqrt_gpu, tau_gpu, np.float32(FREQ), v_gpus[c], grid=grid, block=block, stream=streams[c])
        events[c]['meq'].record(streams[c])
    for c in xrange(CHUNK):
        prev_c = (CHUNK-1+c) % CHUNK
        if c > 0: streams[c].wait_for_event(events[prev_c]['vis'])
        #if c > 0: streams[c].wait_for_event(events[prev_c]['end'])
        # compute vis = dot(v, v.T)
        # transpose below incurs about 20% overhead
        skcuda.cublas.cublasCgemm(h, 't', 'n', NANT, NANT, NPIXC, 1., v_gpus[c].conj().gpudata, NPIXC, v_gpus[c].gpudata, NPIXC, 0., vis_gpus[c].gpudata, NANT)
        events[c]['vis'].record(streams[c])
    for c in xrange(CHUNK):
        prev_c = (CHUNK-1+c) % CHUNK
        if c > 0: streams[c].wait_for_event(events[prev_c]['end'])
        vis_gpus[c].get_async(ary=viss[c], stream=streams[c])
        events[c]['end'].record(streams[c])
    events[CHUNK-1]['end'].synchronize()
    vis = sum(viss)
    print vis
    for c in xrange(CHUNK):
        print c, 'START->END:', events[c]['start'].time_till(events[c]['end']) * 1e-3
        for i,e in enumerate(event_order[:-1]):
            print c, e,'->',event_order[i+1], ':', events[c][e].time_till(events[c][event_order[i+1]]) * 1e-3
    print 'TOTAL:', events[0]['start'].time_till(events[CHUNK-1]['end']) * 1e-3
    #print vis
    ##print np.allclose(aa_cpu, np.dot(a_cpu, a_cpu.T.conj()))

# teardown GPU configuration
skcuda.cublas.cublasDestroy(h)
print 'Done', time.time()

print '\n=== Device attributes'
dev = pycuda.autoinit.device
print 'Name:', dev.name()
print 'Compute capability:', dev.compute_capability()
print 'Concurrent Kernels:', \
     bool(dev.get_attribute(driver.device_attribute.CONCURRENT_KERNELS))
