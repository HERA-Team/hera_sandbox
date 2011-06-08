from time import time
import numpy as np
import aipy as a
import sys,os,optparse

t0 = time()


################################################################
##     Parse inputs
####################################
o = optparse.OptionParser()
o.set_usage('bash_delaycal.py [options] *.ms')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True,src=True)
#o.add_option('--prefix',dest='preprefix',default='PGB',
#    help="""Prepend the output directory with this short thing.  Usually the array
#    acronym [PGB].""")
#o.add_option('--fcat',default='misc',
#    help="Flux calibrator catalog.")
o.add_option('--fconfig',default='120_180_6',
    help='Start_stop_step for output image and cal in MHz.[120_180_6]')
o.add_option('--useflagversion',default=None,type='str',
    help='Version of flags to use [None]')
o.add_option('--minsnr',default=2,type='float',
    help='SNR threshold for calibration solution (see minsnr in CASA bandpass)')
#o.add_option('--scratch',default=None,type='str',
#    help='Directory to use as scratch for imaging. Use to avoid NFS lock errors. Default=None')
#clean off the casa args:
for i in range(len(sys.argv)):
    if sys.argv[i]==inspect.getfile( inspect.currentframe()):break
opts, args = o.parse_args(sys.argv[i+1:])
###################################################


#########################
# user-specified inputs #
#########################
fconfig = n.array(opts.fconfig.split('_')).astype(n.float)
#msfile = '../zen.2455455.67356.uvcbr.ms'
#msfile = vis
#aipycalfile = 'psa455_sim'
aipycalfile = opts.cal#calfile
skymodelcatalogs = opts.cat.split(',')#['misc','three_cr']#['culgoora_160']
ants = '!24'
solint='10min'
uvrange='>20lambda'
fstart  =   fconfig[0]#120
fstop   =   fconfig[1]#180
df      =   fconfig[2]#1
#!0&&1;!2&&3;!4&&5;!6&&7;!8&&9;!10&&11;!12&&13;!14&&15;!16&&17;!18&&19;!20&&21;!22&&23;!24&&25;!26&&27;!28&&29;!30&&31'
docal=True
apply_cal=True
flush = sys.stdout.flush

aa = a.cal.get_aa(aipycalfile,np.array([0.160]))
cat = a.cal.get_catalog(aipycalfile,catalogs = skymodelcatalogs)
cat.compute(aa)

########################
# function definitions #
########################

def aipysrc2dir(aipysrc):
    return me.direction('J2000','%5.6frad'%(aipysrc.ra,),'%5.6frad'%(aipysrc.dec,))

def dir2strlist(d):
    return qa.formxxx(me.getvalue(me.measure(d,'s'))['m0'],format='hms')+' '+qa.formxxx(me.getvalue(me.measure(d,'s'))['m1'],format='dms')

#def how_far(ra1,dec1,ra2,dec2):
#    lambda_diff = ra1 - ra2
#    num = (np.cos(dec2)*np.sin(lambda_diff))**2. + (np.cos(dec1)*np.sin(dec2)-np.sin(dec1)*np.cos(dec2)*np.cos(lambda_diff))**2.
#    denom = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(lambda_diff)
#    return np.arctan2(np.sqrt(num),denom)*(180./3.1415)

def is_src(src,jd,aa):
    return True

#def use_src(cat,Nsrcs=1):
#    visible_cat = {}
#    for src in cat:
#        if is_src(): visible_cat[src] = cat[src]
#    sorted_cat = sorted(visible_cat.iteritems(),key=lambda src: units.pjys(src[1],aa),reverse=True)
#    return ([(src[0],units.pjys(src[1],aa),src[1].ra,src[1].dec,src[1].alt,src[1].index,src[1].mfreq for src in sorted_cat[:Nsrcs]],[dir2strlist(aipysrc2dir(src[1])) for src in sorted_cat[:Nsrcs]])])

def beam(aipysrc,aipycalfile,jd,freq):
    return A
def pjys(src,aa):
    src.compute(aa)
    bm = aa[0].bm_response(src.get_crds('top'))
    return src.jys*bm
def MJDs2JD(t):
    return t/86400 + 2400000.5
def findcalsrcs(aa,cat,Nsrcs=1,altmin=10):
#    ia.open(imagename)
    cat.compute(aa)
    visible_cat = {}
    for src in cat:
#        if srcinimage(ia,cat[src]): visible_cat[src] = cat[src]
        if cat[src].alt>(altmin*n.pi/180):visible_cat[src] = cat[src]
    sorted_cat = sorted(visible_cat.iteritems(),key=lambda src: pjys(src[1],aa),reverse=True)
#    ia.close()
#    return ([(src[0],pjys(src[1],aa),src[1].ra,src[1].dec,src[1].alt,src[1].index,src[1].mfreq) for src in sorted_cat[:Nsrcs]],
#        [dir2strlist(aipysrc2dir(src[1])) for src in sorted_cat[:Nsrcs]])
    return [src[0] for src in sorted_cat[:Nsrcs]]
def dB(A):
    return 10*n.log10(A)
for msfile in args:
    vis=msfile
    ######################
    # naming conventions #
    ######################
    
    cl_name = msfile[:-2]+'cl'
    cal_name = msfile[:-2]+'cal'


    #find the median time
    ms.open(msfile)
    rec = ms.getdata(['time'])
    t = n.median(rec['time'])
    ms.close()
    print "median time: MJD [s], JD [d]"
    print t,MJDs2JD(t)
    aa.set_jultime(MJDs2JD(t))
    
    use_src = ['pic']#,'144']
    print "Choosing a cal source(s)"
    use_src = findcalsrcs(aa,cat,Nsrcs=1,altmin=45)
    print "Using",','.join(use_src)
    
    
    #Make component list
    cl.open()
    for src in use_src:
        s = cat[src]
    #    f = s._jys*(0.16/s.mfreq)**s.index
        f = pjys(cat[src],aa).squeeze()
        w = 5.*(150./160.)
        print '='*50
        print 'Using src:',src
        print 'RA, dec =',s.ra,s.dec
        print 'flux (pjys) =',f
        print 'index =',s.index
        print 'width =',w,'arcmin'
        print '='*50
        flush()
        cl.addcomponent(dir = aipysrc2dir(cat[src]),
                        shape = 'Gaussian',
                        majoraxis = '%2.0farcmin'%w,
                        minoraxis = '%2.0farcmin'%w,
                        positionangle = '0deg',
    #                    flux = [f,0.,0.,0.],
                        flux = f,
                        freq = qa.quantity(160.,'MHz'),
                        spectrumtype = 'spectral index',
                        index = s.index,
                        label = src )
    if os.path.exists(cl_name):
        print "overwriting %s"%cl_name;flush()
        os.system('rm -rf %s'%cl_name)
    cl.rename(cl_name)
    cl.close()
    
    #print "clearcal";flush()
    #clearcal()
    #Make the model
    print "restore known good flags";flush()
    if not opts.useflagversion is None: flagmanager(vis=vis,mode='restore',versionname=opts.useflagversion)
    
    
    print 'Make a Model';flush()
    spw = '0:%d~%dMHz'%(fstart,fstop)
    
    ft(vis=msfile,
        spw=spw,
        complist=cl_name,
        incremental=False)
    
    print '='*50
    pl.figure(10)
    pl.clf()
    pl.figure(11)
    if docal:
        print 'GAINCAL'
        flush()
        bandpass(vis=msfile,
                caltable=cal_name,
                spw=spw,
                selectdata=True,
                timerange='',
                scan='',
                uvrange=uvrange,
                antenna=ants,
                solint=solint,
    #            gaintype='G',
                bandtype='B',
    #            calmode='ap',
                interp=['nearest'],
                fillgaps=10,
                solnorm=False,
                minsnr=opts.minsnr)
        print '='*50
        tb.open(cal_name,nomodify=False)
        G = tb.getcol('GAIN')
        M = tb.getcol('FLAG')
        F = n.linspace(fstart,fstop,num=G.shape[1])
        n.savez(cal_name,G=G[0,:,:],freq=F,mask=M.squeeze())
        lines = []
        #for each spectrum compute a linear delay model
        #then replace the existing channelwise model with the delay model
        dlylog = open(cal_name+'.txt','w')
        for i in range(G.shape[2]):
            P,res,rank,sv,cond = n.ma.polyfit(F/1e3,n.ma.array(
            n.unwrap(n.angle(G[0,:,i]),discont=2.6),mask=M[0,:,i]),1,full=True)
            AP,Ares,Arank,Asv,Acond = n.ma.polyfit(F/1e4,n.ma.array(
            n.abs(G[0,:,i]),mask=M[0,:,i]),2,full=True)
            ampmodel = n.poly1d(AP)
            if rank<2: P,res = [0,0],n.array([0.0])
#            print len(P),P[1],res.squeeze(),rank,cond,sv
            print "Ant: %d,\t Delay [ns]: %3.2f,\t Phase res [r]: %3.2f, \t Amp [Jys/count] %3.2f"%\
            (i,P[1],res.squeeze()/(G.shape[1]-rank),ampmodel(150));flush()
            pl.figure(10)
            l = pl.plot(F,n.ma.masked_where(M[0,:,i],n.unwrap(n.angle(G[0,:,i]),discont=2.6)),label=str(i))[0]
            pl.figure(11)
            pl.plot(F,n.ma.masked_where(M[0,:,i],dB(n.abs(G[0,:,i]/G[0,:,i].max()))),label=str(i),color=l.get_color())
            pl.plot(F,dB(ampmodel(F/1e3)/n.max(ampmodel(F/1e3))),color=l.get_color())
            lines.append(l)
            phasemodel = n.poly1d(P)
            pl.figure(10)
            pl.plot(F,phasemodel(F/1e3),color=l.get_color())
            if apply_cal:
#                G[0,:,i] = n.abs(G[0,:,i])*n.exp(1j*phasemodel(F/1e3))
                G[0,:,i] = ampmodel(F/1e3)*n.exp(1j*phasemodel(F/1e3))
            dlylog.write('%d \t'%i)#output the index
            for p in P:
                dlylog.write('%3.2f\t'%p)
            dlylog.write('%3.2f\t'%(res.squeeze()/(G.shape[1]-rank)))
            dlylog.write('\n')
        dlylog.close()
        if apply_cal:tb.putcol('GAIN',G)
        tb.close()
    #pl.legend(numpoints=1,mode='expand',ncol=8)
    for fi in [10,11]:
        pl.figure(fi)
        ax = pl.gca()
        pl.figlegend(lines,map(str,range(G.shape[2])),'top center',numpoints=1,mode='expand',ncol=8)
        pl.xlabel('Freq [MHz]')
    pl.figure(10)
    pl.ylabel('gain phase [r]')
    pl.savefig(msfile+'.delaymodel.png')
    pl.figure(11)
    pl.ylabel('gain [dB]')
    pl.savefig(msfile+'.ampmodel.png')
    if apply_cal:
        print 'APPLYCAL'
        flush()
        applycal(vis=msfile,
                spw=spw,
    #            selectdata=True,
    #            uvrange=uvrange,
    #            antenna='',
                gaintable=cal_name)
        print '='*50
    print "store the SNR flag table"
    flagmanager(vis=vis,mode='save',versionname='delaycal')
    print 'Computation time: %2.1f m'%((time()-t0)/60.,)
    print '='*50












