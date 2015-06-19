__author__ = 'yunfanzhang'
import optparse, os, sys

#o = optparse.OptionParser()
#o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time point of data')
#o.add_option('-d', '--dft', dest='dst', default=43./3600/24)
#opts,args = o.parse_args(sys.argv[1:])
#print opts, args

def get_fdict(DIR):
    filedict = {}
    args = os.listdir(DIR)
    for s in args:
        s = str(s)
        temp = str(os.path.splittext(s)[0]).split('_')
        t = float(temp[len(temp) - 1])
        filedict[t] = s
    return filedict

def get_file(T, dt, filedict):
    for t in filedict.keys():
        if T-t>0 and T-t < dt: return filedict[t]