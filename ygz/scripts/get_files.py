__author__ = 'yunfanzhang'
import optparse, os, sys, math

#o = optparse.OptionParser()
#o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time point of data')
#o.add_option('-d', '--dft', dest='dst', default=43./3600/24)
#opts,args = o.parse_args(sys.argv[1:])
#print opts, args

#Put files in DIR into dictionary keyed by the fractional julian dates
def get_fdict(DIR):
    filedict = {}
    args = os.listdir(DIR)
    for s in args:
        s = str(s)
        temp = str(os.path.splitext(s)[0]).split('_')
        if temp[0][0] == 'z':             #case of real data zen. simu are pspec_
            lststr = temp[0][4:]                    #strip of zen.
        else: lststr = temp[len(temp) - 1]
        t = float(lststr)
        filedict[t] = s
    return filedict

#Get the filename containing the interval [t,t+dt], containinig T
def get_file(T, dt, filedict):
    for t in filedict.keys():
        if T-t>0 and T-t < dt: return filedict[t]
    print "no file for T=",T, "dt=", dt
    return