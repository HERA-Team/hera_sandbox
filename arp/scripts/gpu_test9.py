#! /home/aparsons/miniconda3/envs/hera_ml/bin/python2.7
from pycuda import compiler, gpuarray, driver
from skcuda.cublas import cublasCreate, cublasSetStream, cublasSgemm, cublasCgemm, cublasDestroy
import numpy as np
from math import ceil
import time
from scipy.interpolate import RectBivariateSpline

NTHREADS = 1024 # make 512 for smaller GPUs
MAX_MEMORY = 2**29 # floats (4B each)

GPU_TEMPLATE = """
// CUDA code for interpolating antenna beams and computing "voltage" visibilities 
// [A^1/2 * I^1/2 * exp(-2*pi*i*freq*dot(a,s)/c)]
// === Template Parameters ===
// "BLOCK_PX": # of sky pixels handled by one GPU block, used to size shared memory
// "NANT"   : # of antennas to pair into visibilities
// "NPIX"   : # of sky pixels to sum over.
// "BEAM_PX": dimension of sampled square beam matrix to be interpolated.  
//            Suggest using odd number.

#include <cuComplex.h>

// Linearly interpolate between [v0,v1] for t=[0,1]
// v = v0 * (1-t) + v1 * t = t*v1 + (-t*v0 + v0)
// Runs on GPU only
__device__
inline float lerp(float v0, float v1, float t) {
    return fma(t, v1, fma(-t, v0, v0));
}

// 3D texture storing beam response on (x=sin th_x, y=sin th_y, nant) grid
// for fast lookup by multiple threads.  Suggest setting first 2 dims of
// bm_tex to an odd number to get pixel centered on zenith.  The pixels
// on the border are centered at [-1,1] respectively.  Note that this
// matrix appears transposed relative to the host-side matrix used to set it.
texture<float, cudaTextureType3D, cudaReadModeElementType> bm_tex;

// Shared memory for storing per-antenna results to be reused among all ants
// for "BLOCK_PX" pixels, avoiding a rush on global memory.
__shared__ float sh_buf[%(BLOCK_PX)s*5];

// Interpolate bm_tex[x,y] at top=(x,y,z) coordinates and store answer in "A"
__global__ void InterpolateBeam(float *top, float *A)
{
    const uint nant = %(NANT)s;
    const uint npix = %(NPIX)s;
    const uint tx = threadIdx.x; // switched to make first dim px
    const uint ty = threadIdx.y; // switched to make second dim ant
    const uint pix = blockIdx.x * blockDim.x + threadIdx.x;
    const uint ant = blockIdx.y * blockDim.y + threadIdx.y;
    const uint beam_px = %(BEAM_PX)s;
    float bm_x, bm_y, px, py, pz, fx, fy, top_z;
    if (pix >= npix || ant >= nant) return;
    if (ty == 0) // buffer top_z for all threads
        sh_buf[tx+%(BLOCK_PX)s * 4] = top[2*npix+pix];
    __syncthreads(); // make sure top_z exists for all threads
    top_z = sh_buf[tx+%(BLOCK_PX)s * 4];
    if (ty == 0 && top_z > 0) { // buffer x interpolation for all threads
        bm_x = (beam_px-1) * (0.5 * top[pix] + 0.5);
        //bm_x = (beam_px) * (0.5 * top[pix] + 0.5) - 0.5f;
        px = floorf(bm_x);   // integer position
        sh_buf[tx+%(BLOCK_PX)s * 0] = bm_x - px; // fx, fractional position
        sh_buf[tx+%(BLOCK_PX)s * 2] = px + 0.5f; // px, pixel index
    }
    if (ty == 1 && top_z > 0) { // buffer y interpolation for all threads
        bm_y = (beam_px-1) * (0.5 * top[npix+pix] + 0.5);
        //bm_y = (beam_px) * (0.5 * top[npix+pix] + 0.5) - 0.5f;
        py = floorf(bm_y);
        sh_buf[tx+%(BLOCK_PX)s * 1] = bm_y - py; // fy, fractional position
        sh_buf[tx+%(BLOCK_PX)s * 3] = py + 0.5f; // py, pixel index
    }
    __syncthreads(); // make sure interpolation exists for all threads
    if (top_z > 0) {
        fx = sh_buf[tx+%(BLOCK_PX)s * 0];
        fy = sh_buf[tx+%(BLOCK_PX)s * 1];
        px = sh_buf[tx+%(BLOCK_PX)s * 2];
        py = sh_buf[tx+%(BLOCK_PX)s * 3];
        pz = ant + 0.5f;
        A[ant*npix+pix] = lerp(lerp(tex3D(bm_tex,px,py,pz),      tex3D(bm_tex,px+1.0f,py,pz),fx),
               lerp(tex3D(bm_tex,px,py+1.0f,pz), tex3D(bm_tex,px+1.0f,py+1.0f,pz),fx), fy);
    } else {
        A[ant*npix+pix] = 0;
    }
    __syncthreads(); // make sure everyone used mem before kicking out
}

// Compute A*I*exp(ij*tau*freq) for all antennas, storing output in v
__global__ void MeasEq(float *A, float *I, float *tau, float freq, cuFloatComplex *v)
{
    const uint nant = %(NANT)s;
    const uint npix = %(NPIX)s;
    const uint tx = threadIdx.x; // switched to make first dim px
    const uint ty = threadIdx.y; // switched to make second dim ant
    const uint row = blockIdx.y * blockDim.y + threadIdx.y; // second thread dim is ant
    const uint pix = blockIdx.x * blockDim.x + threadIdx.x; // first thread dim is px
    float amp, phs;

    if (row >= nant || pix >= npix) return;
    if (ty == 0)
        sh_buf[tx] = I[pix];
    __syncthreads(); // make sure all memory is loaded before computing
    amp = A[row*npix + pix] * sh_buf[tx];
    phs = tau[row*npix + pix] * freq;
    v[row*npix + pix] = make_cuFloatComplex(amp * cos(phs), amp * sin(phs));
    __syncthreads(); // make sure everyone used mem before kicking out
}
"""

def numpy3d_to_array(np_array):
    '''Copy a 3D (d,h,w) numpy array into a 3D pycuda array that can be used 
    to set a texture.  (For some reason, gpuarrays can't be used to do that 
    directly).  A transpose happens implicitly; the CUDA array has dim (w,h,d).'''
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

NTHREADS = 1024 # make 512 for smaller GPUs
MAX_MEMORY = 2**29 # floats (4B each)
MIN_CHUNK = 8

def vis_gpu(antpos, freq, eq2tops, crd_eq, I_sky, bm_cube,
            nthreads=NTHREADS, max_memory=MAX_MEMORY,
            real_dtype=np.float32, complex_dtype=np.complex64,
            verbose=False):
    # ensure shapes
    nant = antpos.shape[0]
    assert(antpos.shape == (nant, 3))
    npix = crd_eq.shape[1]
    assert(crd_eq.shape == (3, npix))
    assert(I_sky.shape == (npix,))
    beam_px = bm_cube.shape[1]
    assert(bm_cube.shape == (nant, beam_px, beam_px))
    ntimes = eq2tops.shape[0]
    assert(eq2tops.shape == (ntimes, 3, 3))
    # ensure data types
    antpos = antpos.astype(real_dtype)
    eq2tops = eq2tops.astype(real_dtype)
    crd_eq = crd_eq.astype(real_dtype)
    Isqrt = np.sqrt(I_sky).astype(real_dtype)
    bm_cube = bm_cube.astype(real_dtype) # XXX complex?
    chunk = max(MIN_CHUNK,2**int(ceil(np.log2(float(nant*npix) / max_memory / 2))))
    npixc = npix / chunk
    # blocks of threads are mapped to (pixels,ants,freqs)
    block = (max(1,nthreads/nant), min(nthreads,nant), 1)
    grid = (int(ceil(npixc/float(block[0]))),int(ceil(nant/float(block[1]))))
    gpu_code = GPU_TEMPLATE % {
            'NANT': nant,
            'NPIX': npixc,
            'BEAM_PX': beam_px,
            'BLOCK_PX': block[0],
    }
    gpu_module = compiler.SourceModule(gpu_code)
    bm_interp = gpu_module.get_function("InterpolateBeam")
    meas_eq = gpu_module.get_function("MeasEq")
    bm_texref = gpu_module.get_texref("bm_tex")
    import pycuda.autoinit
    h = cublasCreate() # handle for managing cublas
    # define GPU buffers and transfer initial values
    bm_texref.set_array(numpy3d_to_array(bm_cube)) # never changes, transpose happens in copy so cuda bm_tex is (BEAM_PX,BEAM_PX,NANT)
    antpos_gpu = gpuarray.to_gpu(antpos) # never changes, set to -2*pi*antpos/c
    Isqrt_gpu = gpuarray.empty(shape=(npixc,), dtype=real_dtype)
    A_gpu = gpuarray.empty(shape=(nant,npixc), dtype=real_dtype) # will be set on GPU by bm_interp
    crd_eq_gpu = gpuarray.empty(shape=(3,npixc), dtype=real_dtype)
    eq2top_gpu = gpuarray.empty(shape=(3,3), dtype=real_dtype) # sent from CPU each time
    crdtop_gpu = gpuarray.empty(shape=(3,npixc), dtype=real_dtype) # will be set on GPU
    tau_gpu = gpuarray.empty(shape=(nant,npixc), dtype=real_dtype) # will be set on GPU
    v_gpu = gpuarray.empty(shape=(nant,npixc), dtype=complex_dtype) # will be set on GPU
    vis_gpus = [gpuarray.empty(shape=(nant,nant), dtype=complex_dtype) for i in xrange(chunk)]
    # output CPU buffers for downloading answers
    vis_cpus = [np.empty(shape=(nant,nant), dtype=complex_dtype) for i in xrange(chunk)]
    streams = [driver.Stream() for i in xrange(chunk)]
    event_order = ('start','upload','eq2top','tau','interpolate','meas_eq','vis','end')
    vis = np.empty((ntimes,nant,nant), dtype=complex_dtype)
    for t in xrange(ntimes):
        if verbose: print '%d/%d' % (t+1, ntimes)
        eq2top_gpu.set(eq2tops[t]) # defines sky orientation for this time step
        events = [{e:driver.Event() for e in event_order} for i in xrange(chunk)]
        for c in xrange(chunk+2):
            cc = c - 1
            ccc = c - 2
            if 0 <= ccc < chunk:
                stream = streams[ccc]
                vis_gpus[ccc].get_async(ary=vis_cpus[ccc], stream=stream)
                events[ccc]['end'].record(stream)
            if 0 <= cc < chunk:
                stream = streams[cc]
                cublasSetStream(h, stream.handle)
                ## compute crdtop = dot(eq2top,crd_eq)
                # cublas arrays are in Fortran order, so P=M*N is actually 
                # peformed as P.T = N.T * M.T
                cublasSgemm(h, 'n', 'n', npixc, 3, 3, 1., crd_eq_gpu.gpudata, 
                    npixc, eq2top_gpu.gpudata, 3, 0., crdtop_gpu.gpudata, npixc)
                events[cc]['eq2top'].record(stream)
                ## compute tau = dot(antpos,crdtop)
                cublasSgemm(h, 'n', 'n', npixc, nant, 3, 1., crdtop_gpu.gpudata, 
                    npixc, antpos_gpu.gpudata, 3, 0., tau_gpu.gpudata, npixc)
                events[cc]['tau'].record(stream)
                ## interpolate bm_tex at specified topocentric coords, store interpolation in A
                ## threads are parallelized across pixel axis
                bm_interp(crdtop_gpu, A_gpu, grid=grid, block=block, stream=stream)
                events[cc]['interpolate'].record(stream)
                # compute v = A * I * exp(1j*tau*freq)
                meas_eq(A_gpu, Isqrt_gpu, tau_gpu, real_dtype(FREQ), v_gpu, 
                    grid=grid, block=block, stream=stream)
                events[cc]['meas_eq'].record(stream)
                # compute vis = dot(v, v.T)
                # transpose below incurs about 20% overhead
                cublasCgemm(h, 'c', 'n', nant, nant, npixc, 1., v_gpu.gpudata, 
                    npixc, v_gpu.gpudata, npixc, 0., vis_gpus[cc].gpudata, nant)
                events[cc]['vis'].record(stream)
            if c < chunk:
                stream = streams[c]
                events[c]['start'].record(stream)
                crd_eq_gpu.set_async(crd_eq[:,c*npixc:(c+1)*npixc], stream=stream)
                Isqrt_gpu.set_async(Isqrt[c*npixc:(c+1)*npixc], stream=stream)
                events[c]['upload'].record(stream)
        events[chunk-1]['end'].synchronize()
        vis[t] = sum(vis_cpus)
        if verbose:
            for c in xrange(chunk):
                print '%d:%d START->END:' % (c, chunk), events[c]['start'].time_till(events[c]['end']) * 1e-3
                #for i,e in enumerate(event_order[:-1]):
                #    print c, e,'->',event_order[i+1], ':', events[c][e].time_till(events[c][event_order[i+1]]) * 1e-3
            print 'TOTAL:', events[0]['start'].time_till(events[chunk-1]['end']) * 1e-3
    # teardown GPU configuration
    cublasDestroy(h)
    return vis

def vis_cpu(antpos, freq, eq2tops, crd_eq, I_sky, bm_cube,
            real_dtype=np.float32, complex_dtype=np.complex64,
            verbose=False):
    nant = len(antpos)
    ntimes = len(eq2tops)
    npix = I_sky.size
    bm_pix = bm_cube.shape[-1]
    Isqrt = np.sqrt(I_sky).astype(real_dtype)
    antpos = antpos.astype(real_dtype)
    A_s = np.empty((nant,npix), dtype=real_dtype)
    vis = np.empty((ntimes,nant,nant), dtype=complex_dtype)
    tau = np.empty((nant,npix), dtype=real_dtype)
    v = np.empty((nant,npix), dtype=complex_dtype)
    bm_pix_x = np.linspace(-1,1,bm_pix)
    bm_pix_y = np.linspace(-1,1,bm_pix)
    for t,eq2top in enumerate(eq2tops.astype(real_dtype)):
        if verbose:
            print '%d/%d' % (t, ntimes)
            t_start = time.time()
        tx,ty,tz = crd_top = np.dot(eq2top, crd_eq)
        for i in xrange(nant):
            spline = RectBivariateSpline(bm_pix_y, bm_pix_x, bm_cube[i], kx=1, ky=1)
            A_s[i] = spline(ty, tx, grid=False)
        A_s = np.where(tz > 0, A_s, 0)
        np.dot(antpos, crd_top, out=tau)
        np.exp(1j*tau, out=v)
        AI_s = A_s * Isqrt
        v *= AI_s
        for i in xrange(len(antpos)):
            np.dot(v[i:i+1].conj(), v[i:].T, out=vis[t,i:i+1,i:])
        if verbose:
            print 'TOTAL:', time.time() - t_start
            #print vis[t].conj()
    np.conj(vis, out=vis)
    return vis

def aa_to_antpos(aa):
    nant = len(aa)
    antpos = np.empty((nant, 3), dtype=np.float32)
    for i in xrange(nant):
        antpos[i] = aa.get_baseline(0,i,src='z')
    return antpos

def aa_to_bm_cube(aa, chan, pol, beam_px=63):
    nant = len(aa)
    assert(pol in ('x','y')) # XXX can only support single linear polarizations right now
    bm_cube = np.empty((nant,beam_px,beam_px), dtype=np.float32) # X is 3rd dim, Y is 2nd dim
    tx = np.linspace(-1,1,beam_px, dtype=np.float32)
    tx = np.resize(tx, (beam_px,beam_px))
    ty = tx.T.copy()
    tx = tx.flatten(); ty = ty.flatten()
    txty_sqr = tx**2 + ty**2
    tz = np.where(txty_sqr < 1, np.sqrt(1-txty_sqr), -1)
    top = np.array([tx,ty,tz])
    aa.select_chans([chan])
    for i in xrange(nant):
        bmi = np.where(tz > 0, aa[i].bm_response(top, pol=pol), 0)
        bm_cube[i] = np.reshape(bmi, (beam_px,beam_px))
    return bm_cube

def aa_to_eq2tops(aa, jds):
    eq2tops = np.empty((len(jds),3,3), dtype=np.float32)
    for i,jd in enumerate(jds):
        aa.set_jultime(jd)
        eq2tops[i] = aa.eq2top_m
    return eq2tops
    
def hmap_to_crd_eq(h):
    px = np.arange(h.npix())
    crd_eq = h.px2crd(px,3).astype(np.float32)
    return crd_eq

def hmap_to_I(h):
    return h[np.arange(h.npix())].astype(np.float32)

if __name__ == '__main__':
    NTIMES = 10
    FREQ = .150

    NANT = 128
    NSIDE = 512
    NPIX = 12*NSIDE**2 # this is assumed to always be bigger than NTHREADS
    CHUNK = max(8,2**int(ceil(np.log2(float(NANT*NPIX) / MAX_MEMORY / 2))))
    NPIXC = NPIX / CHUNK
    BEAM_PX = 63
    GPU = False

    if GPU:
        import time
        import pycuda.autoinit
        print '=== Device attributes'
        dev = pycuda.autoinit.device
        print 'Name:', dev.name()
        print 'Compute capability:', dev.compute_capability()
        print 'Concurrent Kernels:', \
             bool(dev.get_attribute(driver.device_attribute.CONCURRENT_KERNELS))

    print '# Antennas:', NANT
    print 'NSIDE:', NSIDE
    t_start = time.time()
    print 'Starting', t_start
    #print grid, block
    # Initialization of values on CPU side
    np.random.seed(0)
    import aipy
    aa = aipy.cal.get_aa('hsa7458_v001', np.array([FREQ]))
    antpos = aa_to_antpos(aa)
    bm_cube = aa_to_bm_cube(aa, chan=0, pol='x')
    print bm_cube
    #antpos = np.zeros(shape=(NANT,3), dtype=np.float32) # multiply -2pi/c into here
    #antpos[:,0] = 1
    antpos *= -2 * np.pi / 3e1 # cm -> ns conversion
    eq2top = np.array([[1.,0,0],[0,1,0],[0,0,1]])
    eq2tops = np.resize(eq2top, (NTIMES,3,3))
    #crd_eq = np.random.uniform(size=(3,NPIX)).astype(real_dtype)
    crd_eq = np.zeros(shape=(3,NPIX), dtype=np.float32)
    crd_eq[2] = 1
    #crd_eq[0,0] = -1
    crd_eq[0,0] = 1
    # note that bm_tex is transposed relative to the cuda texture buffer
    #bm_cube = np.ones(shape=(NANT,BEAM_PX,BEAM_PX), dtype=np.float32) # X is 3rd dim, Y is 2nd dim
    #bm_cube *= 0
    #bm_cube[:,31,62] = 2
    #bm_cube[:,31,0] = 2
    I_sky = np.zeros(shape=(NPIX,), dtype=np.float32)
    I_sky[0] = 1

    if GPU: vis = vis_gpu(antpos, .15, eq2tops, crd_eq, I_sky, bm_cube, verbose=True)
    else: vis = vis_cpu(antpos, .15, eq2tops, crd_eq, I_sky, bm_cube, verbose=True)
    print np.allclose(vis, 4)
    print 'Time elapsed:', time.time() - t_start
