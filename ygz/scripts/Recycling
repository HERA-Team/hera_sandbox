__author__ = 'yunfanzhang'
# contains recycled code from earlier version

#create a test sample to plot the pairs of points
def test_sample(pairs_final,dt,aa, src,freq,fbmamp,cutoff=5000.): #from select_pair.py
    pairs = []
    bl1,bl2 = pairs_final[0][1][0],pairs_final[0][2][0]
    t1,t2 = pairs_final[0][1][1],pairs_final[0][2][1]
    correlation = pairs_final[0][0]
    print pairs_final[0]
    aa1,aa2 = aa,aa
    src1,src2 = src,src
    while correlation > cutoff:
        aa1.set_jultime(t1)
        src1.compute(aa1)
        alt1 = src1.alt
        aa2.set_jultime(t2)
        src2.compute(aa2)
        alt2 = src2.alt
        if alt1>0 and alt2>0:
            uvw1 = aa1.gen_uvw(*bl1,src=src1).flatten()
            if uvw1[0] < 0: uvw1 = -uvw1
            uvw2 = aa2.gen_uvw(*bl2,src=src1).flatten()
            if uvw2[0] < 0: uvw2 = -uvw2
            pairs.append(((uvw1[0],uvw1[1]),(uvw2[0],uvw2[1])))
            print "sample distance", (uvw1[0]-uvw2[0],uvw1[1]-uvw2[1])
        else: break
        t1,t2 = t1+dt,t2+dt
        correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
    t1,t2 = pairs_final[0][1][1]-dt,pairs_final[0][2][1]-dt
    correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
    while correlation > cutoff:
        aa1.set_jultime(t1)
        src1.compute(aa1)
        alt1 = src1.alt
        aa2.set_jultime(t2)
        src2.compute(aa2)
        alt2 = src2.alt
        if alt1>0 and alt2>0:
            uvw1 = aa1.gen_uvw(*bl1,src=src1).flatten()
            if uvw1[0] < 0: uvw1 = -uvw1
            uvw2 = aa2.gen_uvw(*bl2,src=src1).flatten()
            if uvw2[0] < 0: uvw2 = -uvw2
            pairs.append(((uvw1[0],uvw1[1]),(uvw2[0],uvw2[1])))
            print "sample distance", (uvw1[0]-uvw2[0],uvw1[1]-uvw2[1])
        else: break
        t1,t2 = t1-dt,t2-dt
        correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
    return pairs