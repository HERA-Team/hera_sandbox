import sys, os
import numpy as np
import time
import datetime

class Watch(object):
    def __init__(self,watchFileName='watch.txt',watchoutFile='watch.out'):
        self.watchFileName = watchFileName
        self.watchFileExists = False
        self.watchoutFile = watchoutFile
        self.times = None
        self.ping = None

    def readFile(self,watchFile=None):
        self.watchFileExists = True
        if watchFile is None:
            watchFile = self.watchFileName
        try:
            fp = open(watchFile,'r')
            if self.times is None:
                self.times = []
            if self.ping is None:
                self.ping = []
            elineno = 0
            inEntries = False
            for line in fp:
                line = line.strip()
                if line == 'Watching':
                    elineno = 0
                    inEntries = True
                elif elineno == 1:
                    etime = time.strptime(line,'%a %b %d %H:%M:%S UTC %Y')
                    etime = datetime.datetime.fromtimestamp(time.mktime(etime))
                    #times.append(etime)
                elif elineno == 3:
                    if 'bytes' not in line:
                        pingTime = 0.0
                    else:
                        ds = line.split()
                        try:
                            pingTime = float([x for x in ds if x[0:2]=='ti'][0].split('=')[1])
                        except IndexError:
                            pingTime = 0.0
                    self.ping.append(pingTime)
                    self.times.append(etime)
                if inEntries:
                    elineno+=1
            fp.close()
        except IOError:
            self.watchFileExists = False
            if len(self.times)<1:
                self.times = None
            if len(self.ping)<1:
                self.ping = None

    def writeFile(self,watchoutFile=None):
        if watchoutFile is None:
            watchoutFile = self.watchoutFile
        fp = open(watchoutFile,'w')
        for i,t in enumerate(self.times):
            s = "%s,%.0d\n" % (str(t),self.ping[i])
            fp.write(s)
        fp.close()
        
    def elapsed(self,defaultLong=1.0E6):
        self.diffLast = defaultLong
        if self.times and len(self.times)>0:
            unow = datetime.datetime.utcnow()
            diffLast = unow - self.times[-1]
            self.diffLast = diffLast.total_seconds()/60.0
        return self.diffLast
        

class Alert(object):
    def __init__(self,alertFileName='alerts.txt'):
        self.alertFileExists = False
        self.alertFileName = alertFileName
        self.alerts = None
    def readFile(self,alertFile=None):
        if alertFile is None:
            alertFile = self.alertFileName
        alertFileExists = True
        try:
            fp = open(alertFile,'r')
            if self.alerts is None:
                self.alerts = []
            for line in fp:
                if len(line)>8:
                    line = line.strip().split('.')[0]
                    atime = time.strptime(line,'%Y-%m-%d %H:%M:%S')
                    atime = datetime.datetime.fromtimestamp(time.mktime(atime))
                    self.alerts.append(atime)
            fp.close()
            self.count = np.arange(len(self.alerts))+1
        except IOError:
            self.alertFileExists = False
            if len(self.alerts)<1:
                self.alerts = None

    def elapsed(self,defaultLong=1.0E6):
        self.diffLast = defaultLong
        if self.alerts and len(self.alerts)>0:
            unow = datetime.datetime.utcnow()
            diffLast = unow - self.alerts[-1]
            self.diffLast = diffLast.total_seconds()/60.0
        return self.diffLast
    
    def addEntry(self,unow,alertFile=None):
        if alertFile is None:
            alertFile = self.alertFileName
        fp = open(alertFile,'a')
        fp.write(str(unow)+'\n')
        fp.close()


class Bandwidth(object):
    def __init__(self,measFileName='watchBand.dat',bwFileName='bandwidth.txt',bandTransmitSize=150000000):
        self.measFileName=measFileName
        self.bandTransmitSize=bandTransmitSize
        self.bwFileName = bwFileName

    def getbw(self,bandFile=None):
        if bandFile is None:
            bandFile=self.measFileName
        band = []
        try:
            fp = open(bandFile,'r')
            for line in fp:
               data = line.split()
               band.append(float(data[1]))
        except IOError:
            band = [-1,-1,-1]
        except ValueError:
            band = [-2,-2,-2]
        self.band = band

    def addEntry(self, unow=None, fn='bandwidth.txt'):
        if unow is None:
            unow = datetime.datetime.utcnow()
        fp = open(fn,'a')
        fp.write(str(unow)+' '+str(self.bandTransmitSize)+' ')
        for b in self.band:
            fp.write(str(b)+' ')
        fp.write('\n')
        fp.close()

    def readbwFile(self,bwFile=None):
        if bwFile is None:
            bwFile = self.bwFileName
        bwFileExists = True
        self.times = []
        self.filesize = []
        self.realtime = []
        self.usertime = []
        self.systime = []
        bwData = []
        try:
            fp = open(bwFile,'r')
            for line in fp:
                data = line.split()
                if len(data)==6:
                    lti = data[0]+' '+data[1]
                    lti = lti.strip().split('.')[0]
                    atime = time.strptime(lti,'%Y-%m-%d %H:%M:%S')
                    atime = datetime.datetime.fromtimestamp(time.mktime(atime))
                    self.times.append(atime)
                    self.filesize.append(float(data[2]))
                    self.realtime.append(float(data[3]))
                    self.usertime.append(float(data[4]))
                    self.systime.append(float(data[5]))
            fp.close()
        except IOError:
            self.bwFileExists = False
        self.filesize = np.array(self.filesize)
        self.realtime = np.array(self.realtime)
        self.usertime = np.array(self.usertime)
        self.systime  = np.array(self.systime)

    def makeBigFile(self,bfn=None,size=150,sunit='U'):
        sizes = {'U':1,'kB':1024,'MB':1024*1024,'GB':1024*1024*1024}
        if bfn is None:
            bfn = 'watch%d%s.dat' % (size,sunit)
        with open(bfn, 'wb') as fout:
            fout.write(os.urandom(size*sizes[sunit]))
