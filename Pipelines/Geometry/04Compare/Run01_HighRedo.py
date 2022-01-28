#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'Redo_High'
ID_high = 'High'
ID_redo = 'Redo'
PdbFile = '../../0ConfigData/Pdbs_Under1_good_redo.csv'
PdbDirHigh = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
PdbDirRedo = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
AllGeos = ['N:CA','CA:C','C:N+1','C:O','C-1:N','N:CA:C','C-1:N:CA:C','N:CA:C:N+1','N:N+1','C-1:C','N:O']
GeosHist = ['N:CA','CA:C','C:N+1','C:O','N:CA:C']
Trios = [['C-1:N:CA:C','N:CA:C:N+1','N:CA:C'],['N:CA:C:N+1','N:N+1','N:CA:C'],['N:CA','CA:C','aa'],['C:N+1','C:O','aa']]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import pandas as pd
import A01MakeCsv as A01
import A02CreateHtml as A02
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
make,report = False,True

if make:
    A01.runMakeCsv(ID_high,PdbFile,PdbDirHigh,AllGeos)
    A01.runMakeCsv(ID_redo, PdbFile, PdbDirRedo, AllGeos)

df_geometry_high = pd.read_csv("Csv/" + ID_high + "_01_Geometry.csv")
df_geometry_redo = pd.read_csv("Csv/" + ID_redo + "_01_Geometry.csv")

#Some filtering of the data
aa_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
df_geometry_high = df_geometry_high[df_geometry_high['aa'].isin(aa_list)]
df_geometry_redo = df_geometry_redo[df_geometry_redo['aa'].isin(aa_list)]

df_geometry_high_o = df_geometry_high.query('occupancy == 1')
df_geometry_redo_o = df_geometry_redo.query('occupancy == 1')

df_geometry_high_u = df_geometry_high.query('occupancy < 1')
df_geometry_redo_u = df_geometry_redo.query('occupancy < 1')

if report:
    A02.runCreateReport(ID, ID_high, ID_redo, df_geometry_high_o, df_geometry_redo_o, GeosHist, Trios)
    A02.runCreateReport(ID + '_under', ID_high, ID_redo, df_geometry_high_u, df_geometry_redo_u, GeosHist, Trios)






