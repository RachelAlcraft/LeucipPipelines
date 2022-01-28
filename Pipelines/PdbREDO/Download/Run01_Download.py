
download,convert = True,True
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dataDir = 'C:/Dev/Github/ProteinDataFiles/redo_data/'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
OutFile = '../../0ConfigData/Pdbs_Under1_good_redo.csv'
pdbDir = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
ccp4Dir = "C:/Dev/Github/ProteinDataFiles/ccp4_data_redo/"


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
import A01Download as A01
import A02RunGemmi as A02
import pandas as pd
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


if download:
    A01.download(dataDir, PdbFile, OutFile)

pdb_df = pd.read_csv(OutFile)
pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
pdbs = pdb_df['PDB'].values

if convert:
    A02.runGemmi(pdbs, dataDir, pdbDir, ccp4Dir)


