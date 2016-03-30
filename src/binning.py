import numpy as n

LST_RES = 2*n.pi/24
UV_RES = 1.5

def uv2bin(u,v,lst,uv_res=UV_RES, lst_res=LST_RES):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + int(n.around(lst/lst_res))

def bin2uv(bin, uv_res=UV_RES, lst_res=LST_RES):
    bin = int(bin)
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(lst_res)
    return u,v, lst

def rebin_log(x, y, bin=10):
    '''For y=f(x), bin x into log_e bins, and average y over
    these bin sizes.'''
    logx = n.log(n.abs(x))
    hist1,bins = n.histogram(logx, bins=bin, weights=y)
    hist2,bins = n.histogram(logx, bins=bin)
    logx = .5 * (bins[1:] + bins[:-1])
    return n.e**logx, hist1 / n.where(hist2 == 0, 1., hist2)

def lstbin(lst, lst_res=40.):
    '''Chooses an lst bin for a given lst.  lst_res in seconds'''
    lst_res = lst_res / a.const.sidereal_day * (2*n.pi)
    return bin2uv(uv2bin(0,0,lst,lst_res=lst_res),lst_res=lst_res)[-1]
