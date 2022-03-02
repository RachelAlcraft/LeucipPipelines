#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ID = 'Pairwise'
ID_high = 'PW_High'
ID_redo = 'PW_Redo'
cut_off = 25 # trim outliers
PdbFile = 'C:/Dev/Github/LeucipPipelines/Pipelines/0ConfigData/Pdbs_Under1_good_redo.csv'
PdbDirHigh = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
PdbDirRedo = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
PWAtomsGLY = ['N','CA','C','O','N-1','CA-1','C-1','O-1','N+1','CA+1','C+1','O+1']
PWAtoms = ['CB','CB-1','CB+1']
OtherDssp = ['N:O+2','N:O+3','N:O+4','O:N+2','O:N+3','O:N+4','O:{N}+2','N:{O}+2']
OtherGLY = ['CA-1:C-1:N','C-1:N:CA','N:CA:C','CA:C:N+1','CA:C:O','O:C:N+1','CA-1:CA:CA+1','CA:C:O:N+1']
OtherGLY += ['N+1:CA+1:C+1','N+2:CA+2:C+2','N-1:CA-1:C-1','N-2:CA-2:C-2'] # taus
OtherGLY += ['N:CA:C:N+1','CA-1:C-1:N:CA','CA:C:N+1:CA+1','N:CA:C:O','C-1:N:CA:C','N-1:CA-1:C-1:N']#dihedrals
Other = ['N:CA:CB','CB:CA:C','C-1:N:CA:CB','N:CA:CB:C','CB:CA:C:O','CB:CA:C:N+1']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
aa_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
aa_NO_gly = ['ALA','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import pandas as pd
import sys
from datetime import datetime
sys.path.append('../1Library')
import Helpers as hlp
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def runMakeCsv(ID,PdbFile,PdbDir,AllGeos,aa_filter=[]):
    print('LeucipPipeline: - Creating data -----------------------------------------------------------')
    pdb_df = pd.read_csv(PdbFile)
    pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
    pdbs = pdb_df['PDB'].values
    #pdbs = pdbs[:10]
    ### 2) Loop through and create bio python objects
    print('### Load structures from BioPython #############')
    strucs = bpm.loadPdbStructures(pdbs, PdbDir, extension='ent', prefix='pdb', log=2)

    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(AllGeos, log=1)

    print('### Save dataframe ###############################')
    if len(aa_filter) > 0:
        print('filter on',aa_filter)
        data = data[data['aa'].isin(aa_filter)]

    print('Save to', "Csv/" + ID + "_01_Geometry.csv")
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)


def runMakeCsv_Synthetic(ID,AllGeos):
    import Ext_Geometry as ext_geo
    import Ext_PeptideBuilder as ext_pep
    import A0Class_AtomCollection as rot

    print('Load from', "Csv/" + ID + "_01_Geometry.csv")
    data = pd.read_csv("Csv/" + ID + "_01_Geometry.csv")
    geos = {}
    #dihedrals
    geos['omega'] = 'CA:C:N+1:CA+1'
    geos['phi'] = 'C-1:N:CA:C'
    geos['psi'] = 'N:CA:C:N+1'
    geos['ncaco'] = 'N:CA:C:O'
    geos['cacon1'] = 'CA:C:O:N+1'
    #angles
    geos['tau'] = 'N:CA:C'
    geos['taum1'] = 'C-1:N:CA'
    geos['taup1'] = 'CA:C:N+1'
    geos['caco'] = 'CA:C:O'
    # bonds
    geos['nca'] = 'N:CA'
    geos['cac'] = 'CA:C'
    geos['co'] = 'C:O'
    geos['cn1'] = 'C:N+1'
    # trim the data of all the geos outliers
    print('Length=',len(data.index))
    for geo,atoms in geos.items():#default cut off is 25, so....
        print(geo,min(data[atoms]),max(data[atoms]))
        data = data.sort_values(by=atoms,ascending=True)
        data = data.iloc[cut_off:, :]
        data = data.sort_values(by=atoms, ascending=False)
        data = data.iloc[cut_off:, :]
        print(geo, min(data[atoms]), max(data[atoms]))
    print('Length=', len(data.index))


    print('LeucipPipeline: - Creating synthetic data -----------------------------------------------------------')
    strucs = []
    for i in range(0, 1000):
        structure = ext_pep.initialize_res(geo)
        structure.header = {}
        structure.header['resolution'] = 1
        structure.id = 'X' + str(i)
        for i in range(0, 10):
            #dihedrals
            row = data.sample()
            omega = row[geos['omega']].values[0]
            row = data.sample()
            phi = row[geos['phi']].values[0]
            row = data.sample()
            psi = row[geos['psi']].values[0]
            row = data.sample()
            ncaco =  row[geos['ncaco']].values[0]
            row = data.sample()
            #cacon1 = row[geos['cacon1']].values[0] #library doesn;t allow improer angles
            #angles =
            row = data.sample()
            ncac = row[geos['tau']].values[0]
            row = data.sample()
            taum1 = row[geos['taum1']].values[0]
            row = data.sample()
            taup1 = row[geos['taup1']].values[0]
            row = data.sample()
            caco = row[geos['caco']].values[0]
            # bonds =
            row = data.sample()
            nca = row[geos['nca']].values[0]
            row = data.sample()
            cac = row[geos['cac']].values[0]
            row = data.sample()
            co = row[geos['co']].values[0]
            row = data.sample()
            cn = row[geos['cn1']].values[0]


            print(omega, phi, psi,nca)
            geo = ext_geo.geometry('G')
            #dihedrals =
            geo.omega = float(omega)
            geo.phi = float(phi)
            geo.psi_im1 = float(psi)
            geo.N_CA_C_O_diangle = float(ncaco)

            #angles =
            geo.N_CA_C_angle = float(ncac)
            geo.C_N_CA_angle = float(taum1)
            geo.CA_C_N_angle= float(taup1)
            geo.CA_C_O_angle = float(caco)
            #lengths =
            geo.CA_N_length = float(nca)
            geo.C_O_length = float(co)
            geo.CA_C_length = float(cac)
            geo.peptide_bond = float(cn)

            structure = ext_pep.add_residue(structure, geo)
        strucs.append(structure)

    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(AllGeos, log=1)

    print('### Save dataframe ###############################')
    #if len(aa_filter) > 0:
    #    print('filter on',aa_filter)
    #    data = data[data['aa'].isin(aa_filter)]

    print('Save to', "Csv/PW_" + ID + "_IDEAL_01_Geometry.csv")
    data.to_csv("Csv/" + ID + "_IDEAL_01_Geometry.csv", index=False)



#***************************************************************************************************************
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
##################################################################################
def run(runs):
    start = datetime.now()
    geos = runMakeGeosPairwise(PWAtomsGLY + PWAtoms) + OtherGLY + OtherDssp + Other
    geosGLY = runMakeGeosPairwise(PWAtomsGLY) + OtherGLY + OtherDssp
    for run_for, obs_data, ideal_data in runs:
        print(run_for,obs_data,ideal_data)
        if obs_data:
            if 1 == run_for:
                runMakeCsv(ID_high,PdbFile,PdbDirHigh,geos,aa_filter=aa_NO_gly)
            if 2 == run_for:
                runMakeCsv(ID_high + '_GLY', PdbFile, PdbDirHigh, geosGLY,aa_filter=['GLY'])
            if 3 == run_for:
                runMakeCsv(ID_redo, PdbFile, PdbDirRedo, geos,aa_filter=aa_NO_gly)
            if 4 == run_for:
                runMakeCsv(ID_redo + '_GLY', PdbFile, PdbDirRedo, geosGLY,aa_filter=['GLY'])

        if ideal_data:
            if 2 == run_for:
                runMakeCsv_Synthetic(ID_high + '_GLY', geosGLY)
            if 4 == run_for:
                runMakeCsv_Synthetic(ID_redo + '_GLY', geosGLY)

    end = datetime.now()
    hlp.printTime(start,end)

###########################################################
run([[2,True,True]])
#run([1,False,True],[2,False,True],[3,False,True],[4,False,True])
