#uv = a.spead.UDPUV('192.168.1.1')
uv = a.miriad.UV('file.uv')
ig = make_ig(uv)
for p,d in uv.all():
    ig['spec'] = d
    ...
    ig.sebd
    
    
item group names
#polarization = 'pol'
visibility  = 'vis'
#baseline = 'baseline'
#time = 'time'
flags = 'flags'    
relative baseline coords = 'uvw'
All other items to use miriad variable names
    
