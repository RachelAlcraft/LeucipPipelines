
from LeucipPy import GeoDataFrame as gdf
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Extract information from a data frame
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def getInfo(data, hue,pdb,idcode,id):
    filtered = data.query("pdb_code == '" + pdb + "'")
    idxs = filtered.index
    for idx in idxs:
        thisid = filtered[idcode][idx]
        if thisid.upper() == id.upper():
            return filtered[hue][idx]
    return ''


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# These assist in dataframe manipulation of new keys called y the function, e.g.
# data['atom1'] = data.apply(lambda row: help.applyATOM(row['info' + cfg.GeoA], 0), axis=1)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def getAA(code):
    codes = code.split('|')
    return codes[0],codes[1]

def applyDetailFromOther(detail,other,pos):
    oths = other.split('_')
    if detail == 'GAP':
        a = oths[0]
        b = oths[1]
        aa = a.split('|')
        bb = b.split('|')
        ra = aa[1][:-1]
        rb = bb[1][:-1]
        return abs(int(ra) - int(rb))

    selected_other = oths[pos]
    split_oth = selected_other.split('|')
    aa = split_oth[0]
    rid_chain = split_oth[1]
    atom = split_oth[2]
    if detail == 'RID':
        return rid_chain[:-1]
    if detail == 'CHAIN':
        return rid_chain[-1]
    if detail == 'ATOM':
        return atom
    if detail == 'AA':
        return aa
    return ''
def applyALPHALINK(pdb_code):
    #AF-A0A6M3QH40-F1-model_v2
    afs = pdb_code.split('-')
    if len(afs) > 1:
        af = afs[1]
        lnk = 'https://alphafold.ebi.ac.uk/entry/' + af
        html = '<a target="_blank" href="' + lnk + '">EBI Link</a>'
        return html
    else:
        return pdb_code
def applyGap_rem(other):
    #LYS72A|NZ_HOH661A|O
    oths = other.split('_')
    othsa = oths[0].split('|')
    a,b = getAA(othsa[0])
    ridA = b[:len(b)-1]
    othsb = oths[1].split('|')
    aa, bb = getAA(othsb[0])
    ridB = bb[:len(bb) - 1]
    diff = abs(int(ridA) - int(ridB))
    return diff
def applyRIDDIFF_rem(rid2,rid3):
    return int(rid3) - int(rid2)
def applyRIDGAP_rem(other):
    try:
        oths = other.split('_')
        othso = oths[0].split('|')
        ridO = othso[1][:-1]
        othsa = oths[1].split('|')
        ridA = othsa[1][:-1]
        rO = int(ridO)
        rA = int(ridA)
        ridGap = abs(rO-rA)
        return ridGap
    except:
        return -1
def applyCOORDS(leuDF,atomsdata, pdb_code,info, pos):
    #info = LYS155A|NZ_HOH2127A|O
    #split it into a|b_c|d
    abcd = info.split('_')
    ab = abcd[pos].split('|')
    aa1,x1,atom1 = ab[0],ab[1],ab[2]
    chain1 = x1[-1]
    rid1 = int(x1[:-1])
    #print(info,pdb_code,chain1,aa1,atom1,rid1)
    sp_filtered1 = leuDF.filterDataFrame(atomsdata, inclusions={'pdb_code': [pdb_code], 'chain': [chain1], 'aa': [aa1], 'atom_name': [atom1],'rid': [rid1]})
    if len(sp_filtered1.index) > 0:
        cx, cy, cz = str(round(sp_filtered1['x'].iloc[0],8)), str(round(sp_filtered1['y'].iloc[0],8)), str(round(sp_filtered1['z'].iloc[0],8))
        return cx + ':' + cy + ':' + cz
    else:
        return '::'
def applyCLASS(df,pdb):
    tmpgeo = gdf.GeoDataFrame([])
    filtered = tmpgeo.filterDataFrame(df, inclusions={'PDB': [pdb]})
    vals = filtered['CLASS'].values
    if len(vals) > 0:
        return vals[0]
    return ''
def applyDISTANCE(c1,c2):
    if c1 == '::' or c2 == '::':
        return 0
    from LeucipPy import GeoCalculate as calc
    xyz1 = c1.split(':')
    xyz2 = c2.split(':')
    dis = calc.getDistance(xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2])
    return dis
def applyANGLE(c1,c2,c3):
    if c1 == '::' or c2 == '::' or c3 == '::':
        return 0
    from LeucipPy import GeoCalculate as calc
    xyz1 = c1.split(':')
    xyz2 = c2.split(':')
    xyz3 = c3.split(':')
    dis = calc.getAngle(xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2],xyz3[0],xyz3[1],xyz3[2])
    return dis
def applyID(aa1,rid1,chain1,atom1,aa2,rid2,chain2,atom2,aa3,rid3,chain3,atom3):
    id1 = aa1 + str(rid1) + chain1 + "." + atom1
    id1 += "_" + aa2 + str(rid2) + chain2 + "." + atom2
    id1 += "_" + aa3 + str(rid3) + chain3 + "." + atom3
    return id1

# geranl helper functions
def printTime(start,end,comment=''):
    time_diff = end - start
    print('Time taken (hh:mm:ss.ms) {}'.format(time_diff),comment)
