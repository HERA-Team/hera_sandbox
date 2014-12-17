import datetime, time
import numpy as np

def fmt_time(jd):
    #get the y,m,d out of the jd
    ymd = int(np.floor(jd))
    #from aa.usno.navy.mil
    l = ymd+68569
    n = 4*l/146097
    l = l - (146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l - 1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l
    Y,M,D = i,j,k

    #get h,m,s out of the jd.
    jd -= np.floor(jd)
    jd = 24.*jd
    h = np.floor(jd)
    jd = 60.*(jd-h)
    m = np.floor(jd)
    s = np.floor(60.*(jd-m))

    #format it
    ifmt = '%Y/%m/%d %H:%M:%S'
    tstr = '%d/%d/%d %d:%d:%d'%(Y,M,D,h,m,s)
    t = time.strptime(tstr,ifmt)
    return datetime.datetime(t.tm_year,t.tm_mon,t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec)


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
    
    
    
