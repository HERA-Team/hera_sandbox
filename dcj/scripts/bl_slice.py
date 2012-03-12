import aipy as a,numpy as n,sys
#uv = a.miriad.UV('zen.2455053.97328.uv')
for f in sys.argv[1:]:
#uv = a.miriad.UV('/data1/paper/arp/pgb015/zen.2455020.82248.uv')
    uv = a.miriad.UV(f)
    p,d = uv.read()
    t = p[1]
    uv.rewind()
    data ={}
    for p,d in uv.all():
        if p[1]==t:
            pol = a.miriad.pol2str[uv['pol']]
            bl = a.miriad.bl2ij(uv['baseline'])
            if not data.has_key(pol): data[pol]={bl:n.sum(d[600:800])}
            elif not data[pol].has_key(bl): data[pol][bl]=n.sum(d[600:800])
            else: data[pol][bl].append(n.sum(d[600:800]))
    print n.max(n.array(data['xy'].values())-n.array(data['yx'].values()))
