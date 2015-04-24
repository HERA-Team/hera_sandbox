__author__ = 'yunfanzhang'

#get all baselines equivalent to the given on in antenna array aa

def get_bldict(aa):
    NU,NV = len(aa.ant_layout),len(aa.ant_layout[0])
    nants = len(aa)
    ant_dict, bl_dict = {}, {}
    for i in range(NU):
        for j in range(NV):
            ant_dict[aa.ant_layout[i][j]] = (i,j)  #ant_dict[random ant#]=antlayoutindex
    for i in range(nants):
        for j in range(i+1,nants):
            dkey = (ant_dict[i][0]-ant_dict[j][0], ant_dict[i][1]-ant_dict[j][1])
            if dkey[0]<0 or (dkey[0]==0 and dkey[1]<0): dkey = (-dkey[0],-dkey[1])
            bl_dict[dkey] = bl_dict.get(dkey,[])+[(i,j)]
    return ant_dict, bl_dict

def get_equivalent(ant_dict, bl_dict, ant1, ant2):
    key = (ant_dict[ant1][0]-ant_dict[ant2][0], ant_dict[ant1][1]-ant_dict[ant2][1])
    if key[0]<0 or (key[0]==0 and key[1]<0): key = (-key[0],-key[1])
    return bl_dict[key]

