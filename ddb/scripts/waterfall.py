import numpy as np
import os, sys, os.path
import matplotlib.pyplot as plt
import datetime, time
from jd_cal import *


class wf:
    """Generates waterfall plots using the PAPER flag files
          time1 = datetime of start
          time2 = datetime of stop
          if time1-time2 <= 1 day, plot is in hourse, else it is in days
          reftime = datetime of reference time or 'auto' which makes it equal to the start of the month of time1"""
    
    def __init__(self,time1=datetime.datetime(2012,11,1,0,0,0),time2=datetime.datetime(2012,12,1,0,0,0),reftime='auto',
           Nchan = 1024, band = (100.0,200.0),
           file_path = '/Users/daviddeboer1/Documents/Projects/PAPER/data/rfi/RFI_FILES/2012-2013'):

        #### Get files
        self.ls = []
        for line in os.listdir(file_path):
            dat = line.split('.')
            if dat[-1] == 'npz':
                self.ls.append(line)

        ### Write class variables
        self.time1 = time1
        self.time2 = time2
        self.reftime = reftime
        self.Nchan = Nchan
        self.band = band
        self.file_path = file_path

        #### Set Julian times and frequencies
        self.jd_limits = [get_juldat(time1),get_juldat(time2)]
        if type(reftime) == datetime:
            self.jd_ref = get_juldat(reftime)
        else:
            self.set_autoref()
        self.freqs = np.linspace(band[0], band[1], Nchan)

        self.par()

    def par(self):
        """Prints the run parameters"""
        print '\t--------------------Waterfall--------------------'
        print '\tstart:  %s (%.3f)'% (self.time1.ctime(),self.jd_limits[0])
        print '\tstop:  %s (%.3f)' % (self.time2.ctime(),self.jd_limits[1])
        if type(self.reftime) == datetime:
            print '\tref:  %s (%.3f)'  % (self.reftime.ctime(),self.jd_ref)
        else:
            print '\tref:  auto  at  %s (%.3f)' % (self.reftime,self.jd_ref)
        print '\tband:  %.3f - %.3f' % (self.band[0],self.band[1])
        print '\tchannels:  %d' % (self.Nchan)


    def time(self,lim='start', year=2012, month=11, day=1, hour=0, minute=0, second=0):
        """Set the start, stop and reference times"""
        if lim[0].lower()=='r' and type(year) == str:
            self.set_autoref()
            t = self.reftime
            jd = self.jd_ref
        else:
            if year<2000:
                year+=2000
            t = datetime.datetime(year,month,day,hour,minute,second)
            jd = get_juldat(t)
            
        if lim[0:3].lower() == 'sta' or lim[0].lower() == 'b':
            self.time1 = t
            self.jd_limits[0] = jd
            self.set_autoref()
            lim = 'start'
        elif lim[0:3].lower() == 'sto' or lim[0].lower() == 'e':
            self.time2 = t
            self.jd_limits[1] = jd
            lim = 'stop'
        elif lim[0].lower() == 'r':
            self.reftime = t
            self.jd_ref = jd
            lim = 'reference'
        else:
            print lim+' time not understood'
            return 0
        print 'Setting %s time to %s (%f)' % (lim,str(t),jd)
        return t, jd

    def set_autoref(self):
        """sets reference time in 'auto' based on time1 (start)"""
        t = datetime.datetime(self.time1.year,self.time1.month,1,0,0,0)
        self.reftime = t.ctime()
        self.jd_ref = get_juldat(t)
        print 'Updated autoref to '+self.reftime

    def run(self):
        """script-like run"""
        self.get_wf()
        self.timePlot()
        self.scanPlot()
        
        
    def get_wf(self,periodName=None):
        """Reads the files within the time window and computes the waterfall.
                periodName is the name of the time period - if none uses start day as 'Mon_D_YYYY'
                Data stored in:
                    self.wftimes - time of the scanline
                    self.wfplot - fraction by channel within a file
                    self.aveFlags - fraction by channel over whole time period (gets plotted and written out)
           scanPlot() plots the taken data line-by-line (with lines for long breaks)
           timePlot() plots them as a function of civil time from the reference time"""
        #### Initialize variables
        flags = None
        aveFlags = None
        times = []   
        wfplot = []    ### this is the waterfall plot image
        wftimes = []   ### ...and these are the times
        check = {}     ###-- this is used to keep track of Nchans, just in case

        ###--------------------Loop over files to get data
        firstOne = True
        stopNow = False
        periodName = None
        num_pts_total = 0
        for fial in self.ls:
            flags = None
            fial_day = float(fial.split('.')[1] + '.' + fial.split('.')[2])
            if fial_day-self.jd_limits[0] < -0.1:
                continue
            dat = np.load(os.path.join(self.file_path,fial))
            filetime = 0.0  # used to get "average" time of file
            num_pts_in_file = 0
            ###------------Loop over times in file
            for i in range(len(dat['times'])):
                dattime = dat['times'][i]
                ###--- keep track of Nchan
                flagkey = len(dat[str(i)])
                if check.has_key(len(dat[str(i)])):
                    check[flagkey]+=1
                else:
                    check[flagkey]=0
                ###---
                ###-----Check if right number of channels and in time window
                if len(dat[str(i)]) == self.Nchan and dattime >= self.jd_limits[0] and dattime <= self.jd_limits[1]:
                    if firstOne:
                        print 'Starting in file '+fial
                        asctime = datetime.datetime.ctime(fmt_time(dattime))
                        firstOne = False
                        if periodName is None:
                            periodName = asctime.split()[1]+'_'+asctime.split()[2]+'_'+asctime.split()[-1]
                        print asctime
                    times.append(dattime)
                    filetime+=dattime
                    num_pts_in_file+=1
                    if flags == None:
                        flags = dat[str(i)].astype(np.float)
                    else:
                        flags += dat[str(i)]
                    if aveFlags == None:
                        aveFlags = dat[str(i)].astype(np.float)
                    else:
                        aveFlags += dat[str(i)]
                elif dattime >= self.jd_limits[1]:
                    print 'Stopping at file '+fial
                    print datetime.datetime.ctime(fmt_time(dattime))
                    stopNow = True
                    break
                ###-----
            if flags is not None:
                num_pts_total += num_pts_in_file
                flags /= num_pts_in_file
                filetime /= num_pts_in_file
                wftimes.append(filetime)
                wfplot.append(flags)
                plt.plot(self.freqs,flags)
            dat.close()
            if stopNow:
                break
            ###----------
        aveFlags /= num_pts_total
        outfile = periodName+'_aveFlag.txt'
        fp = open(outfile,'w')
        for i in range(len(self.freqs)):
            s = '%f\t%f\n' % (self.freqs[i],aveFlags[i])
            fp.write(s)
        fp.close()
        ###--------------------
        plt.figure('average')
        plt.plot(self.freqs,aveFlags)

        tnorm = self.jd_ref
        self.times = np.array(times) - tnorm
        self.wftimes = np.array(wftimes) - tnorm
        self.wfplot = np.array(wfplot)
        self.aveFlags = aveFlags
        self.periodName = periodName
        return periodName

    
    #-------------------- "grid" data
    def evenTimes(self,dial=2.0,upper=15.0,justNights=False,plot=False,clr=1.0,verbose=True):
        """Puts the waterfall plot on a more or less even time grid, adding null scanlines as necessary (set to clr)
           The new waterfall data is now stored in:
               self.t
               self.wf
            Some diagnostic plots and verbosity etc"""
        wftimes = self.wftimes
        wfplot = self.wfplot
        
        d2m = 24.0*60.0
        if verbose:
            print 'Eventimes dial: %f, upper:  %f' % (dial,upper)
        interval = 0.0
        nint = 0
        for k in range(len(wftimes)-1):
            pint = (wftimes[k+1] - wftimes[k])
            if pint < upper/d2m:
                interval+=pint
                nint+=1
        interval/=(nint)
        print 'interval is ',interval*d2m
        t = []
        wf = []
        addedInterval = 0
        for i in range(len(wftimes)-1):
            checkInterval = abs((wftimes[i+1] - wftimes[i])/interval)
            t.append(wftimes[i])
            wf.append(wfplot[i])
            if checkInterval >= dial:  # we've skipped some times...so need to add in those blanks as clr
                addedInterval+=1
                if justNights:
                    nAdd = 4
                else:
                    nAdd = int(np.floor( checkInterval ) ) - 1
                    if nAdd == 0:
                        nAdd = 1
                if not justNights and verbose:
                    print '\taddedInterval %d between %f to %f' % (addedInterval,wftimes[i]*d2m,wftimes[i+1]*d2m)
                    print '\t\tcheckInterval',checkInterval
                    print '\t\tnAdd = ',nAdd
                for j in range(nAdd):
                    newEntry = wftimes[i]+(j+1)*interval
                    t.append(newEntry)
                    wf.append(np.zeros(np.shape(wfplot[0])) + clr)
        t.append(wftimes[-1])
        wf.append(wfplot[-1])
        t = np.array(t)
        wf = np.array(wf)
        dint = []
        for i in range(len(t)-1):
            dt = t[i+1] - t[i]
            dint.append(dt)
        dint.append(interval)
        dint = np.array(dint)
        maxInterval = dint.max()
        minInterval = dint.min()
        if verbose:
            print 'maxInterval = ',maxInterval*d2m
            print 'minInterval = ',minInterval*d2m
        if plot:
            plt.figure(plot)
            plt.plot(t,dint*d2m)
        self.t = np.array(t)
        self.wf = np.array(wf)
        return (minInterval,interval,maxInterval)

    def timePlot(self,dial=1.6,upper=15.0,save=True):
        """plot "timed-based" waterfall and save=True/False.  Does the iterative work to put waterfall on even timeline"""

        t_original = np.copy(self.wftimes)
        p_original = np.copy(self.wfplot)

        interval = self.evenTimes(dial=dial,upper=upper,plot=False,verbose=False)
        maxInterval = interval[2]
        aveInterval = interval[1]
        print 'Starting intervals:  %.4f, %.4f, %.4f' % (interval[0],interval[1],interval[2])
        while maxInterval > aveInterval*dial:
            interval = self.evenTimes(dial=dial,upper=upper,plot=1000,verbose=False)
            maxInterval = interval[2]
            aveInterval = interval[1]
            self.wftimes = np.copy(self.t)
            self.wfplot = np.copy(self.wf)
        print 'Ending intervals:  %.4f, %.4f, %.4f' % (interval[0],interval[1],interval[2])
        self.wftimes = np.copy(t_original)
        self.wfplot = np.copy(p_original)
        
        if (self.wftimes[-1] - self.wftimes[0]) <= 1.0:
            xl = 'Hour of day '+str(np.floor(self.t[0]))
            self.t = (self.t - np.floor(self.t[0]))*24.0
        else:
            xl = 'Day of month'
        plt.figure(2000)
        plt.title(self.periodName)
        plt.imshow(self.wf,origin='lower',extent=[100.0,200.0,self.t[0],self.t[-1]])
        plt.axis('normal')
        plt.axis([self.freq[0],self.freq[-1],np.floor(self.t[0]),np.ceil(self.t[-1])])
        plt.xlabel('Frequency [MHz]')
        plt.ylabel(xl)
        plt.colorbar()
        plt.show()
        if save:
            plt.savefig(self.periodName+'_time.png',dpi=300)

    def scanPlot(self,dial=10.0,upper=15.0,save=True):
        """ plot "scan-based" waterfall and save=True/False.  Puts in gaps for long time gaps (ie day-time)"""
        interval = self.evenTimes(dial=dial,upper=upper,justNights=True,plot=False,verbose=False)
        plt.figure(2001)
        plt.title(self.periodName)
        plt.imshow(self.wf,origin='lower',extent=[100.0,200.0,0,len(self.t)])
        plt.axis('normal')
        plt.axis([self.freq[0],self.freq[-1],0,len(self.t)])
        plt.xlabel('Frequency [MHz]')
        plt.colorbar()
        plt.show()
        if save:
            plt.savefig(self.periodName+'_scan.png',dpi=300)


    def otherPlotsEtc(self):
        """plot some other stuff"""
        plt.figure(2)
        plt.plot(self.times,self.times,'.')
        plt.plot(self.wftimes,self.wftimes,'ro')
        plt.plot(t,t,'gx')



