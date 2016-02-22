import aipy
import capo
import numpy

def sep2mdlbl():
    ''' Returns a dictionary indexed by separation type (ex: "0,2"), containing the corresponding model baselines that Omnical outputs in its model visibilities'''
    npz1 = numpy.load('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_xtalk/zen.2456942.38963.xx.uvcRRE.npz')
    mdl = {}
    #keys are sep ("0,2") where the first number is how many antennas to move down the array and second number is how many antennas to move right 
    #values are vismdl keys ("<1,4> xx")
    aa = aipy.cal.get_aa('psa6622_v003',numpy.array([.15]))
    sep, conj = capo.red.group_redundant_bls(aa.ant_layout)
    for key in npz1.keys():
        if key[0] == "<":
            i,j = key.split(',')[0].split('<')[1],key.split(',')[1].split(' ')[0].split('>')[0]
            bl = aipy.miriad.ij2bl(int(i),int(j))
            for s in sep.values(): #s is array of similar baselines
                if bl in s:
                    sepkey = sep.keys()[sep.values().index(s)]
            mdl[sepkey] = key.split(' ')[0].split(',')[0].split('<')[1]+'_'+key.split(' ')[0].split(',')[1].split('>')[0]
    return mdl

def mdlbl2sep():
    ''' Returns a dictionary indexed by model baseline (ex: "1_4"), containing the corresponding separation type ("0,2")'''
    mdl = sep2mdlbl()
    seps = {}
    for sep in mdl.keys():
        seps[mdl[sep]] = sep
    return seps
