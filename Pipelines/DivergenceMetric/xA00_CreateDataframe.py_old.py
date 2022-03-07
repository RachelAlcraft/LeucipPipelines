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
import A0Class_AtomCoordinates as coords
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
    # trim the data of all the geos outliers
    for pair_rama in [0,1,2]:
        csv_file = "Csv/" + ID + "_IDEAL_UNPAIRED_01_Geometry.csv"
        if pair_rama == 1:
            csv_file = "Csv/" + ID + "_IDEAL_PAIRED_01_Geometry.csv"
        elif pair_rama == 2:
            csv_file = "Csv/" + ID + "_IDEAL_SAMPLED_01_Geometry.csv"

        print('LeucipPipeline: - Creating synthetic data',csv_file,' -----------------------------------------------------------')
        strucs = []
        for i in range(0, 1000):
            structure = ext_pep.initialize_res(ext_geo.geometry('G'))
            structure.header = {}
            structure.header['resolution'] = 1
            structure.id = 'X' + str(i)
            last_psi = None
            for i in range(0, 10):
                param_gen = coords.AtomCoordinates(data, log=1)
                if pair_rama == 1:
                    geom, psi_next = param_gen.generateRamaGeo()
                elif pair_rama == 2:
                    geom,psi_next = param_gen.generateSampleGeo()
                else:
                    geom, psi_next = param_gen.generateAllRandomGeo()

                if last_psi != None:
                    geom.psi_im1 = last_psi
                structure = ext_pep.add_residue(structure, geom,last_psi)
            strucs.append(structure)

        print('### Creating dataframe for correlations #############')
        geo_mak = dfm.GeometryMaker(strucs, log=1)
        data = geo_mak.calculateGeometry(AllGeos, log=1)

        print('### Save dataframe ###############################')
        #if len(aa_filter) > 0:
        #    print('filter on',aa_filter)
        #    data = data[data['aa'].isin(aa_filter)]

        print('Save to', csv_file)
        data.to_csv(csv_file, index=False)



#***************************************************************************************************************

##################################################################################
def run(runs):
    import A0_Globals as globals
    start = datetime.now()
    geos = globals.getGeos(False)
    geosGLY = globals.getGeos(True)
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
run([[2,False,True]])
#run([1,False,True],[2,False,True],[3,False,True],[4,False,True])
