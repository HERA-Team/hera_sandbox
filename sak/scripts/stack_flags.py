import numpy as np
import pylab
import glob
from matplotlib.ticker import NullFormatter
nullfmt   = NullFormatter()

flagdata = np.zeros((3912, 1024))
percent_flagged_freq = np.zeros((1024))
percent_flagged_int = np.zeros((3912))

c=0.
for npzFile in glob.glob('*_rfidata.npz'):
	dat = np.load(npzFile)
	
	tempflags = dat['flags']
	print 'shape of "flag waterfall" ',npzFile,tempflags.shape
	
	temp_pcnt_freq = dat['pcnt_freq']
	temp_pcnt_int = dat['pcnt_int']
	
	if tempflags.shape!=flagdata.shape: continue
	c+=1.
	flagdata += tempflags
	percent_flagged_freq += temp_pcnt_freq
	percent_flagged_int += temp_pcnt_int
	
print 'Num files:',c


####Plotting
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02
rect_imshow = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

pylab.figure(1, figsize=(8,8))		

axImshow = pylab.axes(rect_imshow)
axHistx = pylab.axes(rect_histx)
axHisty = pylab.axes(rect_histy)
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

axImshow.imshow(flagdata/c,aspect='auto',interpolation='nearest') #AVERAGING -- NOTE THE "/C"

axHistx.plot(range(percent_flagged_freq.shape[0]),percent_flagged_freq/c)#AVERAGING -- NOTE THE "/C"
axHisty.plot(percent_flagged_int/c,range(percent_flagged_int.shape[0]))#<< & vv AVERAGING -- NOTE THE "/C" 

axHistx.fill_between(range(percent_flagged_freq.shape[0]),0,percent_flagged_freq/c,color='blue',alpha=0.3)

axHistx.set_xlim( axImshow.get_xlim() )
axHisty.set_ylim( axImshow.get_ylim() )
axHisty.set_xlim( (0.,100.) )

axHistx.set_ylabel(r'% $\rm{flagged}$ $\nu$')
axHisty.set_xlabel(r'% $\rm{flagged}$ $t$')

axImshow.set_xlabel(r'Frequency bin')
axImshow.set_ylabel(r'Integration')

pylab.show()
