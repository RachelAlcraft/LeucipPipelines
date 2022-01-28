#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
make,filter, randomise, html_boundaries = False,False,False,True
run_for = [1] #High, High_GLY, REDO,. REDO_GLY


ID = 'Pairwise'
ID_high = 'PW_High'
ID_redo = 'PW_Redo'
PdbFile = '../../0ConfigData/Pdbs_Under1_good_redo.csv'
PdbDirHigh = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
PdbDirRedo = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
PWAtomsGLY = ['N','CA','C','O','N-1','CA-1','C-1','O-1','N+1','CA+1','C+1','O+1']
PWAtoms = ['CB','CB-1','CB+1']
OtherGLY = ['O:(N,O,S@2)','O:(N,O,S@3)','O:(N,O,S@4)','O:(N,O,S@5)','O:(N,O,S@6)','O:(N,O,S@7)','O:{N}+2','N:{O}+2','C-1:N:CA','N:CA:C','CA:C:N+1','CA:C:O','O:C:N+1','N:CA:C:N+1','CA-1:CA:CA+1','CA:C:O:N+1','C-1:N:CA:C','CA-1:C-1:N:CA','CA:C:N+1:CA+1']
Other = ['N:CA:CB','CB:CA:C','C-1:N:CA:CB','N:CA:CB:C','CB:CA:C:O','CB:CA:C:N+1']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
geosHist = ['N:CA','CA:C','C:O','C:N+1','N:CA:C','C-1:N:CA','CA:C:N+1','N:CA:C:N+1']
#triples = [['N:CA:C:N+1','N:N+1','N:CA:C',100,120],['C-1:N:CA:C','C:C-1','N:CA:C',100,120],['CA:C:N+1','N:CA:C','C:O',None,None],['C:O','O:{N}+2','CA:C:N+1',None,None],['C:O','O:(N,O,S)+2','HB_atom',None,None],['C:O','O:(N,O,S)+2','HB_aa',None,None],['O:(N,O,S)+2','O:(N,O,S@2)+2','C:O',None,None],['O:(N,O,S)+2','O:(N,O,S@2)+2','HB_atom',None,None],['O:(N,O,S)+2','O:(N,O,S@2)+2','HB_atom2',None,None]]
triples = [['N:CA:C:N+1','N:N+1','N:CA:C',100,120]]
triples.append(['N:CA:C:N+1','N:N+1','PROB',100,120])
triples.append(['C-1:N:CA:C','C:C-1','N:CA:C',100,120])
triples.append(['C-1:N:CA:C','C:C-1','PROB',100,120])
triples.append(['CA:C:N+1','N:CA:C','C:O',None,None])
triples.append(['CA:C:N+1','O:CB','C:O',None,None])
triples.append(['CA:C:N+1','C:CB','C:O',None,None])
triples.append(['O:CB','C:O','CA:C:N+1',None,None])
triples.append(['C:CB','C:O','CA:C:N+1',None,None])
triples.append(['C:CB','C:O','N:CA:C',None,None])
triples.append(['N:CA:CB:C','C:CB','N:CA:C',None,None])
triples.append(['N:CA:CB:C','N:CA:C','CA:C:N+1',None,None])
triples.append(['N:CA:CB:C','N:CA:C','@CA:C:N+1',None,None])
triples.append(['N:CA:CB:C','N:CA:C','C:O',None,None])
triples.append(['N:CA:CB:C','CA:C:N+1','C:O',None,None])
triples.append(['N:CA:CB:C','CA:C:N+1','N:C-1',None,None])
triples.append(['N:CA:CB:C','CB:CA:C','N:CA:C',None,None])
triples.append(['N:CA:CB:C','N:CA:C','CA+1:CB',None,None])
triples.append(['N:CA:CB:C','N:CA:C','N:CA:C:N+1',None,None])
triples.append(['N:CA:CB:C','N:CA:C','PROB',None,None])
triples.append(['N:CA:CB:C','N:CA:C','@N:CA:C:N+1',None,None])

triples.append(['N:CA:C','N:CA:CB','CB:CA:C',None,None])
triples.append(['CB:CA:C','N:CA:C','N:CA:CB',None,None])
triples.append(['N:CA:CB','CB:CA:C','N:CA:C',None,None])

triples.append(['N:CA:C','N:CA:CB','CA:C:N+1',None,None])
triples.append(['CB:CA:C','N:CA:C','CA:C:N+1',None,None])
triples.append(['N:CA:CB','CB:CA:C','CA:C:N+1',None,None])

triples.append(['N:CA:C','N:CA:CB','C:O',None,None])
triples.append(['CB:CA:C','N:CA:C','C:O',None,None])
triples.append(['N:CA:CB','CB:CA:C','C:O',None,None])
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
#Some filtering of the data
aa_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
aa_NO_gly = ['ALA','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

psi_boundaries = [-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]
tau_boundaries = [110,112,114,116]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import pandas as pd
import sys
from datetime import datetime
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Helpers as hlp
import A01MakeCsv as A01
import A02CreateHtml as A02
import A03FilterAndRandomise as A03

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
start = datetime.now()
sets= ['PW_High', 'PW_High_GLY','PW_Redo', 'PW_Redo_GLY']

if make:
    geos = A01.runMakeGeosPairwise(PWAtomsGLY + PWAtoms) + OtherGLY + Other
    geosGLY = A01.runMakeGeosPairwise(PWAtomsGLY) + OtherGLY
    if 1 in run_for:
        A01.runMakeCsv(ID_high,PdbFile,PdbDirHigh,geos,aa_filter=aa_NO_gly)
    if 2 in run_for:
        A01.runMakeCsv(ID_high + '_GLY', PdbFile, PdbDirHigh, geosGLY,aa_filter=['GLY'])
    if 3 in run_for:
        A01.runMakeCsv(ID_redo, PdbFile, PdbDirRedo, geos,aa_filter=aa_NO_gly)
    if 4 in run_for:
        A01.runMakeCsv(ID_redo + '_GLY', PdbFile, PdbDirRedo, geosGLY,aa_filter=['GLY'])

if filter:
    print('Loading original dataframes and filtering')
    if 1 in run_for:
        df_geometry_high = pd.read_csv("Csv/" + ID_high + "_01_Geometry.csv")
        df_geometry_high = df_geometry_high.query('occupancy == 1')
        df_geometry_high = df_geometry_high.query('bfactor <= 10')
        df_geometry_high['HB_atom1'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@2)'], 1), axis=1)
        df_geometry_high['HB_aa1'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('AA', row['infoO:(N,O,S@2)'], 1), axis=1)
        df_geometry_high['HB_atom2'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@3)'], 1), axis=1)
        df_geometry_high['HB_atom3'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@4)'], 1), axis=1)
        df_geometry_high['HB_atom4'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@5)'], 1), axis=1)
        df_geometry_high['HB_atom5'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@6)'], 1), axis=1)
        df_geometry_high['HB_atom6'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@7)'], 1), axis=1)
        df_geometry_high.to_csv("Csv/" + ID_high + "_02_Geometry.csv", index=False)
        A03.filterDataframes('PSI', 'N:CA:C:N+1','PW_High', df_geometry_high, psi_boundaries)
        A03.filterDataframes('TAU','N:CA:C','PW_High', df_geometry_high, tau_boundaries)
    if 2 in run_for:
        df_geometry_high_gly = pd.read_csv("Csv/" + ID_high + "_GLY_01_Geometry.csv")
        df_geometry_high_gly = df_geometry_high_gly.query('occupancy == 1')
        df_geometry_high_gly = df_geometry_high_gly.query('bfactor <= 10')
        df_geometry_high_gly['HB_atom1'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@2)'], 1), axis=1)
        df_geometry_high_gly['HB_aa1'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('AA', row['infoO:(N,O,S@2)'], 1), axis=1)
        df_geometry_high_gly['HB_atom2'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@3)'], 1), axis=1)
        df_geometry_high_gly['HB_atom3'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@4)'], 1), axis=1)
        df_geometry_high_gly['HB_atom4'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@5)'], 1), axis=1)
        df_geometry_high_gly['HB_atom5'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@6)'], 1), axis=1)
        df_geometry_high_gly['HB_atom6'] = df_geometry_high_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@7)'], 1), axis=1)
        df_geometry_high_gly.to_csv("Csv/" + ID_high + '_GLY' + "_02_Geometry.csv", index=False)
        A03.filterDataframes('PSI', 'N:CA:C:N+1', 'PW_High_GLY', df_geometry_high_gly, psi_boundaries)
        A03.filterDataframes('TAU', 'N:CA:C', 'PW_High_GLY', df_geometry_high_gly, tau_boundaries)
    if 3 in run_for:
        df_geometry_redo = pd.read_csv("Csv/" + ID_redo + "_01_Geometry.csv")
        df_geometry_redo = df_geometry_redo.query('occupancy == 1')
        df_geometry_redo = df_geometry_redo.query('bfactor <= 10')
        df_geometry_redo['HB_atom'] = df_geometry_redo.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S)+2'], 1), axis=1)
        df_geometry_redo['HB_aa'] = df_geometry_redo.apply(lambda row: hlp.applyDetailFromOther('AA', row['infoO:(N,O,S)+2'], 1), axis=1)
        df_geometry_redo['HB_atom2'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@2)+2'], 1), axis=1)
        df_geometry_redo.to_csv("Csv/" + ID_redo + "_02_Geometry.csv", index=False)
        A03.filterDataframes('PSI', 'N:CA:C:N+1', 'PW_Redo', df_geometry_redo, psi_boundaries)
        A03.filterDataframes('TAU', 'N:CA:C', 'PW_Redo', df_geometry_redo, tau_boundaries)
    if 4 in run_for:
        df_geometry_redo_gly = pd.read_csv("Csv/" + ID_redo + "_GLY_01_Geometry.csv")
        geometry_redo_gly = df_geometry_redo_gly.query('occupancy == 1')
        df_geometry_redo_gly = df_geometry_redo_gly.query('bfactor <= 10')
        df_geometry_redo_gly['HB_atom'] = df_geometry_redo_gly.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S)+2'], 1), axis=1)
        df_geometry_redo_gly['HB_aa'] = df_geometry_redo_gly.apply(lambda row: hlp.applyDetailFromOther('AA', row['infoO:(N,O,S)+2'], 1), axis=1)
        df_geometry_redo_gly['HB_atom2'] = df_geometry_high.apply(lambda row: hlp.applyDetailFromOther('ATOM', row['infoO:(N,O,S@2)+2'], 1), axis=1)
        df_geometry_redo_gly.to_csv("Csv/" + ID_redo + '_GLY' + "_02_Geometry.csv", index=False)
        A03.filterDataframes('PSI', 'N:CA:C:N+1', 'PW_Redo_GLY', df_geometry_redo_gly, psi_boundaries)
        A03.filterDataframes('TAU', 'N:CA:C', 'PW_Redo_GLY', df_geometry_redo_gly, tau_boundaries)

if randomise:
    #randomise PSI set
    for i in range(0,4):
        set_nm = sets[i]
        if i+1 in run_for:
            set_id = 'PSI'
            bound_num = len(psi_boundaries)
            A03.randomise(set_id, set_nm, bound_num)

            # randomise TAU set
            set_id = 'TAU'
            bound_num = len(tau_boundaries)
            A03.randomise(set_id, set_nm, bound_num)

if html_boundaries:
    for i in range(0, 4):
        st = sets[i]
        if i+1 in run_for:
            #original
            print(st)
            df_high = pd.read_csv("Csv/" + st + "_02_Geometry.csv")
            #psi randomised
            df_high_psi_rand = pd.read_csv("Csv/PSI_" + st + "_randomised_04_Geometry.csv")
            #tau randomised
            df_high_tau_rand = pd.read_csv("Csv/TAU_" + st + "_randomised_04_Geometry.csv")

            df_list = [[st,df_high],[st + '_PSI_RAND',df_high_psi_rand],[st + '_TAU_RAND',df_high_tau_rand]]

            A02.runCreateBoundariesReport(st +'_Bound', df_list, geosHist, triples)

end = datetime.now()
hlp.printTime(start,end)


