import aipy as a, numpy as n, sys
try: import atpy
except(ImportError): pass
import numpy as np
import re
"""
Codes from danny
"""
def pol2ind(pol):
    assert(pol<-4)
    return pol+8

def make_indexable(x):
    try: len(x)
    except(TypeError):x.shape = (1,)
    return x
def dB(x):
    return 10*n.log10(x)
def idB(x):
    return 10**(x/10)
def select_table_where_source(table,srcname):
    #Assumes MRC names!!!
    #make a special exception for pic
    if srcname=='pic':
        seq=20658
    elif srcname=='cen':
        seq=61767
    else:
        result = table.where(table['Name']== 'MRC %s'%srcname)
        seq = result['Seq']
    return table.where(table['Seq']==seq)
def spectrum(table,fmax=10e3):
    #return freqs,fluxes,errors below fmax (default 10e3 MHz)
    freqs,fluxes,errors = n.array(table['nu']),n.array(table['S_nu_']/1e3),n.array(table['e_S_nu_']/1e3)
    freqs = make_indexable(make_indexable(freqs)[freqs<fmax])
    fluxes = make_indexable(make_indexable(fluxes)[freqs<fmax])
    errors = make_indexable(make_indexable(errors)[freqs<fmax])
    return freqs,fluxes,errors
def ned_spectrum(nedxmlfile,doprint=False,fmax=5000e6):
    ############
    ### set some data from NED
    ####
    NED = atpy.Table(nedxmlfile)
    fluxcol = 'photo_col7'
    freqcol = 'photo_col6'
    errorcol = 'photo_col8'
    unitscol = 'photo_col5'
    refcol = 'photo_col10'
    NEDfreqs = NED[freqcol]
    NED = NED.where(NEDfreqs<fmax)
    nedfreqs = n.array(NED[freqcol]/1e6)
    nedfluxes = n.array(NED[fluxcol])
    nederrorstrings = NED[errorcol]
    nedrefs = NED[refcol]
    #nederrors = n.array([0.]*len(nedfluxes)) %Danny being conservative. No error bar = somehow suspect measurement.
    #    in this past has indicated duplicates data using an older calibration
    nederrors = nedfluxes*0.25 #arp being liberal, assuming 25% error bars for data with no catalog error
    guess = n.ones_like(nedfluxes).astype(n.bool)
    for i,err in enumerate(nederrorstrings):
        if len(err)<2: continue
        try:
            nederrors[i] = float(err.split('+/-')[1])
            guess[i] = False
        except(ValueError): 
            if err.split('+/-')[1].endswith('%'):
                nederrors[i] = nedfluxes[i]*float(err.split('+/-')[1].split('%')[0].strip())/100
                guess[i] = False
    if doprint:
        print "using %s with the following error bars"%nedxmlfile
        print "freq\tflux\tError\tFrac\tGuess?\tCite"
        for i in range(len(nederrors)):
            print "%7.2f\t%7.2f\t%7.2f\t%7.2f\t%r\t%s"%(
                    nedfreqs[nederrors>0][i],
                    nedfluxes[nederrors>0][i],
                    nederrors[nederrors>0][i],
                    nederrors[nederrors>0][i]/nedfluxes[nederrors>0][i],
                    guess[nederrors>0][i],
                    nedrefs[nederrors>0][i])
                    

    return nedfreqs[nederrors>0],nedfluxes[nederrors>0],nederrors[nederrors>0]
def find_trace_peak(trace,gridsize):
    #input trace and number of grid points
    #output indices and values
    H,bins = np.histogramdd(trace)
    peak_index = H.argmax()
    return peak_index,[bins[i][peak] for i,j in enumerate(peak_index)]
def find_closest(A,a):
    return np.abs(A-a).argmin()
def find_percentile_errors(trace,percentile,nbins=100):
    thist,bins = np.histogram(trace,bins=nbins)
    binwidth = np.diff(bins)[0]
    lower = bins[find_closest(np.cumsum(thist)/float(np.sum(thist)),
        0.5-(percentile/100.)/2)]+binwidth/2.
    upper = bins[find_closest(np.cumsum(thist)/float(np.sum(thist)),
        0.5+(percentile/100.)/2)]+binwidth/2.
    med = np.median(trace)
    return [lower,med,upper]
def get_votable_column_names(tablefile,tid=0):
    from atpy.votable import parse
    votable = parse(tablefile)
    colnames = {}
    for id, table in enumerate(votable.iter_tables()):
        if id==tid:
            break    
    for field in table.fields:
        colnames[field._ID] = field.name
    return colnames
def load_table(tablefile):
    colnames = get_votable_column_names(tablefile)
    table = atpy.Table(tablefile)
    for ID,name in colnames.iteritems():
        try:
            table.rename_column(ID,name)
        except(Exception):
            continue
    return table
def a2l(array):
    return ','.join(map(str,array))
def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])

def condition_goodness(X,goodness=.9):
    #checks the condition of a matrix
    #return true if its reasonably invertible
    #return false if the condition approaches the limit
    # set by the bit dept of the input array
    #
    #the 'goodness' parameter is defined as the fractional number of
    # bits in X's precision above which the matrix is assumed to be
    #ill conditioned
    # based on: http://mathworld.wolfram.com/ConditionNumber.html
    condition = np.linalg.cond(X)
    bit_count = X.itemsize*8
    if np.iscomplexobj(X):
        bit_count /= 2 
    bit_goodness_threshold = np.round(bit_count * goodness)
    return np.log(condition)/np.log(2)  < bit_goodness_threshold

    

