import numpy as np

def count_flags(m, labels):
    order = np.argsort(np.sum(m, axis=1))[::-1]
    m = m[order][:,order]
    labels = [labels[i] for i in order]
    cnts = {}
    for i,label in enumerate(labels): cnts[label] = np.sum(m[i,i:])
    return cnts

def check_ants(reds, data, flag_thresh=.2):
    all_bls = reduce(lambda x,y: x+y, reds)
    auto_pwr, ant2col = {}, {}
    for bl in all_bls:
        try: d = data[bl]
        except(KeyError): d = data[bl[::-1]].conj()
        auto_pwr[bl] = np.median(np.sum(np.abs(d)**2, axis=0))
        ant2col[bl[0]] = ant2col.get(bl[0],len(ant2col))
        ant2col[bl[1]] = ant2col.get(bl[1],len(ant2col))
    col2ant = {}
    for i in ant2col: col2ant[ant2col[i]] = i
    nants = len(ant2col)
    Fant = np.zeros((nants,nants), dtype=np.int)
    for bls in reds:
        C = np.zeros((len(bls),len(bls)), dtype=np.float)
        for i,bl1 in enumerate(bls):
            try: d1 = data[bl1]
            except(KeyError): d1 = data[bl1[::-1]].conj()
            for j,bl2 in enumerate(bls[i:]):
                j += i
                try: d2 = data[bl2]  
                except(KeyError): d2 = data[bl2[::-1]].conj()
                pwr12 = np.median(np.abs(np.sum(d1*d2.conj(), axis=0)))
                C[i,j] = C[j,i] = pwr12 / np.sqrt(auto_pwr[bl1] * auto_pwr[bl2])
        cnts = count_flags(np.where(C < flag_thresh, 1, 0), bls)
        for i,j in cnts: Fant[ant2col[i],ant2col[j]] = Fant[ant2col[j],ant2col[i]] = cnts[(i,j)]
    cnts = count_flags(np.where(Fant >= 1, 1, 0), [col2ant[i] for i in xrange(nants)])
    return cnts 
