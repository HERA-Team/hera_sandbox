def convert(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    return [map(float, option.split('_')) for option in arg]

def gen_chans(chanopt, uv, coords, is_delay):
    """Return an array of active channels and whether or not a range of
    channels is selected (as opposed to one or more individual channels)
    based on command-line arguments."""
    global import numpy as n
    is_chan_range = True
    if chanopt == 'all': chans = n.arange(uv['nchan'])
    else:
        chanopt = convert(chanopt)
        if coords != 'index':
            if is_delay:
                def conv(c):
                    return int(n.round(c * uv['sdf'] * uv['nchan'])) \
                        + uv['nchan']/2
            else:
                def conv(c): return int(n.round((c - uv['sfreq']) / uv['sdf']))
        else:
            if is_delay:
                def conv(c): return int(c) + uv['nchan']/2
            else:
                def conv(c): return c
        chanopt = [map(conv, c) for c in chanopt]
        if len(chanopt[0]) != 1: 
            chanopt = [n.arange(x,y, dtype=n.int) for x,y in chanopt]
        else: is_chan_range = False
        chans = n.concatenate(chanopt)
    return chans.astype(n.int), is_chan_range
