#input a 3d casa image cube and bin it by df
#imagename = 
#df = #in Hz

h = imhead(imagename=imagename,mode='list')
freqs = n.arange(0,h['shape'][2])*h['cdelt3']+h['crval3']
freq_bins = n.arange(freqs.min(),freqs.max(),df)
axes='freq'
function='mean'
for i in range(len(freq_bins)-1):
    #print "freq range:",i,freq_bins[i],freq_bins[i+1]
    mfreq = (freq_bins[i+1] + freq_bins[i])/2/1e6
    outfile=imagename.replace('sm2arcmin','sm2arcmin_%dMHz'%mfreq)
    print outfile
    chanmin = n.abs(freqs-freq_bins[i]).argmin()
    chanmax = n.abs(freqs-freq_bins[i+1]).argmin()
    chans='%d~%d'%(chanmin,chanmax)
    imcollapse()
    

    

