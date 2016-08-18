import aipy as a, numpy as n

def read_files(filenames, antstr, polstr, decimate=1, decphs=0, verbose=False, recast_as_array=True):
    info = {'lsts':[], 'times':[]}
    ts = {}
    dat, flg = {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        if decimate > 1: uv.select('decimate', decimate, decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if not ts.has_key(t):
                info['times'].append(t)
                info['lsts'].append(uv['lst'])
                ts[t] = None
            bl = (i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            pol = a.miriad.pol2str[uv['pol']]
            if not dat[bl].has_key(pol):
                dat[bl][pol],flg[bl][pol] = [],[]
            dat[bl][pol].append(d)
            flg[bl][pol].append(f)
    info['freqs'] = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    if recast_as_array:
        # This option helps reduce memory footprint, but it shouldn't
        # be necessary: the replace below should free RAM as quickly
        # as it is allocated.  Unfortunately, it doesn't seem to...
        for bl in dat.keys():
          for pol in dat[bl].keys():
            dat[bl][pol] = n.array(dat[bl][pol])
            flg[bl][pol] = n.array(flg[bl][pol])
        info['lsts'] = n.array(info['lsts'])
        info['times'] = n.array(info['times'])
    return info, dat, flg

