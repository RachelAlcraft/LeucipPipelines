
aa_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
aa_GLY = {'aa':['GLY']}
aa_NO_gly = {'aa':['ALA','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']}

def runMakeGeosPairwise(atoms_list):
    #Make a pairwise listout of atoms
    geos = []
    for i in range(0, len(atoms_list)):
        atmA = atoms_list[i]
        for j in range(i+1, len(atoms_list)):
            atmB = atoms_list[j]
            if atmA != atmB:
                geos.append(atmA + ':' + atmB)
    return geos

def getGeos(onlyGLY):
    PWAtomsGLY = ['N', 'CA', 'C', 'O', 'N-1', 'CA-1', 'C-1', 'O-1', 'N+1', 'CA+1', 'C+1', 'O+1']
    PWAtoms = ['CB', 'CB-1', 'CB+1']
    OtherDssp = ['N:O+2', 'N:O+3', 'N:O+4', 'O:N+2', 'O:N+3', 'O:N+4', 'O:{N}+2', 'N:{O}+2']
    OtherGLY = ['CA-1:C-1:N', 'C-1:N:CA', 'N:CA:C', 'CA:C:N+1', 'CA:C:O', 'O:C:N+1', 'CA-1:CA:CA+1', 'CA:C:O:N+1']
    OtherGLY += ['N+1:CA+1:C+1', 'N+2:CA+2:C+2', 'N-1:CA-1:C-1', 'N-2:CA-2:C-2']  # taus
    OtherGLY += ['N:CA:C:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1', 'N:CA:C:O', 'C-1:N:CA:C','N-1:CA-1:C-1:N','N:CA-1:C-1:O-1']  # dihedrals
    Other = ['N:CA:CB', 'CB:CA:C', 'C-1:N:CA:CB', 'N:CA:CB:C', 'CB:CA:C:O', 'CB:CA:C:N+1']

    if onlyGLY:
        return runMakeGeosPairwise(PWAtomsGLY) + OtherGLY + OtherDssp
    else:#not using CB for now
        #return runMakeGeosPairwise(PWAtomsGLY + PWAtoms) + OtherGLY + OtherDssp + Other
        return runMakeGeosPairwise(PWAtomsGLY) + OtherGLY + OtherDssp

def getGeosToAbs():
    return ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1','N:CA-1:C-1:O-1']

def excludeGeos():
    return ['N-2:CA-2:C-2','N-3:CA-3:C-3']

def trimData(data, trim,geos):
    # only stabdard amino acids
    aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    temp_data = data[data['aa'].isin(aa_list)]
    temp_data = temp_data[temp_data['aa-1'].isin(aa_list)]
    temp_data = temp_data[temp_data['aa+1'].isin(aa_list)]
    for geo in geos:
        temp_data = temp_data.sort_values(by=geo,ascending=True)
        temp_data = temp_data.iloc[trim:,:]
        temp_data = temp_data.sort_values(by=geo, ascending=False)
        temp_data = temp_data.iloc[trim:, :]
    geos_to_abs = getGeosToAbs()
    for gabs in geos_to_abs:
        temp_data[gabs] = abs(temp_data[gabs])
    temp_data = temp_data.query('`CA:C:N+1:CA+1` >= 100')
    temp_data = temp_data.query('`CA-1:C-1:N:CA` >= 100')

    return temp_data

