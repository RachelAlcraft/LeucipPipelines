
from LeucipPy import GeoDataFrame as gdf
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Extract information from a data frame
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def getInfo(data, hue,pdb,id):
    filtered = data.query("pdb_code == '" + pdb + "'")
    idxs = filtered.index
    for idx in idxs:
        thisid = filtered['ID'][idx]
        if thisid.upper() == id.upper():
            return filtered[hue][idx]
    return ''


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# These assist in dataframe manipulation of new keys called y the function, e.g.
# data['atom1'] = data.apply(lambda row: help.applyATOM(row['info' + cfg.GeoA], 0), axis=1)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def applyOther(other):
    oths = other.split('_')
    otha = oths[1][:3]
    othb =oths[1].split('|')[1]
    return otha + "_" + othb
def applyGap(other):
    #LYS72A|NZ_HOH661A|O
    oths = other.split('_')
    othsa = oths[0].split('|')
    ridA = othsa[0][3:len(othsa[0])-1]
    othsb = oths[1].split('|')
    ridB = othsb[0][3:len(othsb[0]) - 1]
    diff = abs(int(ridA) - int(ridB))
    return diff
def applyAA(other):
    #LYS72A|NZ_HOH661A|O
    oths = other.split('_')
    otha = oths[1][:3]
    return otha
def applyRIDDIFF(rid2,rid3):
    return int(rid3) - int(rid2)
def applyCHAIN(other):
    oths = other.split('_')
    othsa = oths[1].split('|')
    ch = othsa[0][len(othsa[0]) - 1:]
    return str(ch)
def applyRID(other):
    oths = other.split('_')
    othsa = oths[1].split('|')
    ridA = othsa[0][3:len(othsa[0]) - 1]
    return ridA
def applyATOM(other,pos):
    oths = other.split('_')
    othsa = oths[pos].split('|')
    atA = othsa[1]
    return atA
def applyCOORDS(leuDF,atomsdata, pdb_code,info, pos):
    #info = LYS155A|NZ_HOH2127A|O
    infos = info.split('_')
    info1 = infos[pos]
    info1s = info1.split('|')
    chain1 = info1s[0][len(info1s[0]) - 1:]
    aa1 =  info1[:2]
    if info1[2] in ['0','1','2','3','4','5','6','7','8','9']:
        #if aa1 not in ['ZN','MG','CO','CU','NA','CA','FE']:
        rid1 = int(info1s[0][2:len(info1s[0]) - 1])
        aa1 = info1[:2]
    else:
        aa1 = info1[:3]
        rid1 = int(info1s[0][3:len(info1s[0]) - 1])



    atom1 = info1s[1]
    sp_filtered1 = leuDF.filterDataFrame(atomsdata, inclusions={'pdb_code': [pdb_code], 'chain': [chain1], 'aa': [aa1], 'atom_name': [atom1],'rid': [rid1]})
    if len(sp_filtered1.index) > 0:
        cx, cy, cz = str(round(sp_filtered1['x'].iloc[0],4)), str(round(sp_filtered1['y'].iloc[0],4)), str(round(sp_filtered1['z'].iloc[0],4))
        return cx + ':' + cy + ':' + cz
    else:
        print(info)
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
def applyID(aa1,rid1,chain1,atom1,aa2,rid2,chain2,atom2,aa3,rid3,chain3,atom3):
    id1 = aa1 + str(rid1) + chain1 + "." + atom1
    id1 += "_" + aa2 + str(rid2) + chain2 + "." + atom2
    id1 += "_" + aa3 + str(rid3) + chain3 + "." + atom3
    return id1

