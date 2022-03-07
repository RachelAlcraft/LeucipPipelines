#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ID = 'Pairwise'
ID_high = 'PW_High'
ID_redo = 'PW_Redo'
cut_off = 25 # trim outliers
PdbFile = 'C:/Dev/Github/LeucipPipelines/Pipelines/0ConfigData/Pdbs_Under1_good_redo.csv'
PdbDirHigh = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
PdbDirRedo = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
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
def runMakeCsv(ID,PdbFile,PdbDir,AllGeos,aa_filter):
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
    print('filter on', aa_filter)
    if aa_filter != {}:
        print('filter on',aa_filter,len(data.index))
        data = geo_mak.filterDataFrame(data, {'aa': ['GLY']})
        print('filterred', len(data.index))


    print('Save to', "Csv/" + ID + "_01_Geometry.csv")
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)

#***************************************************************************************************************

##################################################################################
def run(runs):
    import A0_Globals as globals
    start = datetime.now()
    geos = globals.getGeos(False)
    geosGLY = globals.getGeos(True)
    for run_for in runs:
        print(run_for)
        if 1 == run_for:
            runMakeCsv(ID_high,PdbFile,PdbDirHigh,geos,aa_filter=globals.aa_NO_gly)
        if 2 == run_for:
            runMakeCsv(ID_high + '_GLY', PdbFile, PdbDirHigh, geosGLY,aa_filter=globals.aa_GLY)
        if 3 == run_for:
            runMakeCsv(ID_redo, PdbFile, PdbDirRedo, geos,aa_filter=globals.aa_NO_gly)
        if 4 == run_for:
            runMakeCsv(ID_redo + '_GLY', PdbFile, PdbDirRedo, geosGLY,aa_filter=globals.aa_GLY)
    end = datetime.now()
    hlp.printTime(start,end)

###########################################################
run([2])
#run([1,False,True],[2,False,True],[3,False,True],[4,False,True])
