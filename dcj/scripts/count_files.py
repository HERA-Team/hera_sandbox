path='/data3/PAPER/psa128/'
from glob import glob
import numpy as n
import os
NIGHTS = glob(path+'/psa*')
count = {}
totalcount = 0
start_PJD = 6620
times_PJD = []
for night in NIGHTS:
    FILES = glob(night+'/*RRE')
    count[os.path.basename(night)]=len(FILES)
    totalcount += len(FILES)
    try:
        night_PJD = int(night[-4:])
        times_PJD.append(night_PJD)
    except:
        continue
times_PJD = n.array(times_PJD)
FILES_EXPECTED = (times_PJD.max() - times_PJD.min())*280
print "total files found at ASU",totalcount
print "fraction of season ",float(totalcount)/(FILES_EXPECTED)*100
