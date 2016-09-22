import numpy as np
import matplotlib.pyplot as plt
import os
import datetime

# Note that after 2456874.62465  the 'labjack' temperature got added
added_labjack_into_list = 2456874.62465  #don't use this for the moment
# but I do add 1 to all tmap entries from 'bE' and up
tmap = {'jd':0,'labjack':1,'bE':2,'cW':3,'bW':4,'cE':5,'r1a':9,'r1b':10,'r2a':11,'r2b':12,
        'r3a':13,'r3b':14,'r4a':16,'r4b':17,'r5a':18,'r5b':19,'r6a':20,'r6b':21,
        'r7a':23,'r7b':24,'r8a':25,'r8b':26}
K2C = 273.15
alla = ['r1a','r2a','r3a','r4a','r5a','r6a','r7a','r8a']
allb = ['r1b','r2b','r3b','r4b','r5b','r6b','r7b','r8b']
ae = ['bE','cE']
aw = ['bW','cW']
r1 = ['r1a','r1b']
r2 = ['r2a','r2b']
r3 = ['r3a','r3b']
r4 = ['r4a','r4b']
r5 = ['r5a','r5b']
r6 = ['r6a','r6b']
r7 = ['r7a','r7b']
r8 = ['r8a','r8b']

def pt(fn,Temps,pltnum=0):
    d = np.loadtxt(fn)
    plt.figure(pltnum)
    for t in Temps:
        plt.plot(d[:,tmap['jd']],d[:,tmap[t]]-K2C,label=t)
    plt.legend()

def all(Temps,pltnum=0,dayspast=10.0,legend=False):
    if type(Temps) == str:
        Temps = list(Temps)
    original_datalen = 29
    now = datetime.datetime.utcnow()
    print now
    reftime = datetime.datetime(now.year,now.month,now.day,0,0,0)
    refjd = get_juldat(reftime)
    
    Files = []
    for line in os.listdir('.'):
        if line[0:6] == 'temp.2':
            data = line.split('.')
            acc = pow(10.0,len(data[2]))
            fjd = float(data[1]) + float(data[2])/acc
            if refjd - fjd < dayspast:
                Files.append(line)
    Files.sort()
    temp = []
    ll1 = 0
    print 'Using '+Files[0]+' - '+Files[-1]
    for f in Files:
        fs = os.stat(f).st_size
        if fs<100000:
            print 'File too small:  %s  (%d)' % (f,fs)
            continue
        #print 'Opening %s  (%d bytes)' % (f,fs)

        fp = open(f,'r')
        for line in fp:
            data = line.split()
            ddd = []
            for d in data:
                ddd.append(float(d))
            if len(ddd) != ll1:
                print 'Change in len to %d in %s  (%d)' % (len(ddd),f,fs)
                ll1 = len(ddd)
            if len(ddd) == original_datalen:  #add in a 0 value for the non-existent labjack term
                ddd[1:1] = [0.0]
            #    for i in range(default_datalen - len(ddd)):
            #        ddd.append(0.0)
            temp.append(ddd)
    temp = np.array(temp)

    temp[:,0] = 24.0*(temp[:,0]-refjd)
    plt.figure(pltnum)
    for t in Temps:
        col = tmap[t]
        plt.plot(temp[:,tmap['jd']],temp[:,col]-K2C,'.',label=t)
    if legend:
        plt.legend(loc='lower left')
    #plt.axis(xmin=0.0,xmax=24.0)
    s = 'Hours since %d-%d-%d (utc)' % (now.year,now.month,now.day)
    plt.title(s)
    plt.xlabel('Time [h]')
    plt.ylabel('Degrees [C]')
    
    return temp


def get_caldat(jd):
    """From Astronomy on a Personal Computer, Montenbruck and Pfleger"""
    a = long(jd+0.5)
    if a < 2299161:   ### Julian calendar
        b = 0
        c = a + 1524
    else:               ### Gregorian
        b = long((a-1867216.25)/365.25)
        c = a + b - (b/4) + 1525
    d = long( (c-122.1)/365.25 )
    e = 365*d + d/4
    f = long( (c-e)/30.6001 )
    Day = c - e - int(30.6001*f)
    Month = f - 1 - 12*(f/14)
    Year = d - 4715 - ((7+Month)/10)
    FracOfDay = jd+0.5 - np.floor(jd+0.5)
    Hour = 24.0*FracOfDay
    Minute = 60.0*(Hour - np.floor(Hour))
    Hour = np.floor(Hour)
    Second = 60.0*(Minute - np.floor(Minute))
    Minute = np.floor(Minute)
    Microsecond = 1000.0*(Second - np.floor(Second))
    Second = np.floor(Second)
    return datetime.datetime(Year,Month,Day,Hour,Minute,Second,Microsecond)

def get_juldat(time):
    """From Astronomy on a Personal Computer, Montenbruck and Pfleger"""
    Year = time.year
    Month = time.month
    Day = time.day
    Hour = time.hour
    Minute = time.minute
    Second = time.second
    if Month <= 2:
        Month+=12
        Year-=1
    if 10000*Year + 100*Month + Day <= 15821004:
        b = -2 + ((Year+4716)/4) - 1179             ### Julian calendar
    else:
        b = (Year/400) - (Year/100) + (Year/4)      ### Gregorian calendar
    MjdMidnight = 365*Year - 679004 + b + int(30.6001*(Month+1)) + Day
    FracOfDay = (Hour + Minute/60.0 + Second/3600.0)/24.0

    jd = 2400000.5 + MjdMidnight + FracOfDay

    return jd
