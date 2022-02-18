#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ID = 'Pairwise'
ID_high = 'PW_High'
ID_redo = 'PW_Redo'
PdbFile = 'C:/Dev/Github/LeucipPipelines/Pipelines/0ConfigData/Pdbs_Under1_good_redo.csv'
PdbDirHigh = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
PdbDirRedo = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
PWAtomsGLY = ['N','CA','C','O','N-1','CA-1','C-1','O-1','N+1','CA+1','C+1','O+1']
PWAtoms = ['CB','CB-1','CB+1']
OtherDssp = ['N:O+2','N:O+3','N:O+4','O:N+2','O:N+3','O:N+4','O:{N}+2','N:{O}+2']
OtherGLY = ['C-1:N:CA','N:CA:C','CA:C:N+1','CA:C:O','O:C:N+1','N:CA:C:N+1','CA-1:CA:CA+1','CA:C:O:N+1','C-1:N:CA:C','CA-1:C-1:N:CA','CA:C:N+1:CA+1']
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


def runMakeCsv_Synthetic(ID,AllGeos,aa_filter=[]):
    import Ext_Geometry as ext_geo
    import Ext_PeptideBuilder as ext_pep
    import A0Class_AtomCollection as rot

    print('LeucipPipeline: - Creating synthetic data -----------------------------------------------------------')
    omega_rules = rot.RotationRules('0.5{-180,-150}:0.5{150,180}')
    phi_rules = rot.RotationRules('0.5{-180,-50}:0.5{50,180}')
    psi_rules = rot.RotationRules('4{-180,-120}:1{-120,-50}:4{-50,50}:1{50,120}:4{120,180}')
    strucs = []
    for i in range(0, 1000):
        omega = int(omega_rules.getRandomRotation())
        phi = int(phi_rules.getRandomRotation())
        psi = int(psi_rules.getRandomRotation())
        #print(omega, phi, psi)
        geo = ext_geo.geometry('G')
        geo.omega = omega
        geo.phi = phi
        geo.psi_im1 = psi
        structure = ext_pep.initialize_res(geo)
        structure.header = {}
        structure.header['resolution'] = 1
        structure.id = 'X' + str(i)

        for i in range(0, 10):
            structure = ext_pep.add_residue(structure, geo)
        strucs.append(structure)

    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(AllGeos, log=1)

    print('### Save dataframe ###############################')
    #if len(aa_filter) > 0:
    #    print('filter on',aa_filter)
    #    data = data[data['aa'].isin(aa_filter)]

    print('Save to', "Csv/PW_" + ID + "_01_Geometry.csv")
    data.to_csv("Csv/PW_" + ID + "_01_Geometry.csv", index=False)



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
def run(run_for):
    start = datetime.now()
    geos = runMakeGeosPairwise(PWAtomsGLY + PWAtoms) + OtherGLY + OtherDssp + Other
    geosGLY = runMakeGeosPairwise(PWAtomsGLY) + OtherGLY + OtherDssp
    if 1 in run_for:
        runMakeCsv(ID_high,PdbFile,PdbDirHigh,geos,aa_filter=aa_NO_gly)
    if 2 in run_for:
        runMakeCsv(ID_high + '_GLY', PdbFile, PdbDirHigh, geosGLY,aa_filter=['GLY'])
    if 3 in run_for:
        runMakeCsv(ID_redo, PdbFile, PdbDirRedo, geos,aa_filter=aa_NO_gly)
    if 4 in run_for:
        runMakeCsv(ID_redo + '_GLY', PdbFile, PdbDirRedo, geosGLY,aa_filter=['GLY'])
    if 5 in run_for:
        runMakeCsv_Synthetic('SYN_GLY', geosGLY, aa_filter=['G'])

    end = datetime.now()
    hlp.printTime(start,end)

###########################################################
run([5])
