import pandas as pd

unzip, csv, html  = False, False, True
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
geos = ['N:N+1', 'CA:CA+1','C:C+1','N:O','O:CA+1']
pairs = [['N:N+1', 'CA:CA+1'],['N:N+1','C:C+1'],['N:N+1','N:O'],['N:N+1','O:CA+1'],['CA:CA+1','C:C+1'],['CA:CA+1','N:O'],['CA:CA+1','O:CA+1'],['C:C+1','N:O'],['C:C+1','O:CA+1'],['N:O','O:CA+1']]
chunk = 200
name='OneFour'
################################################
pdbDir = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import A01Unzip as A01
import A02CreateCsv as A02
import A03CreateHtml as A03

pdb_df = pd.read_csv(PdbFile)
pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
pdbs = pdb_df['PDB'].values


if unzip:
    A01.unzipProteome(pdbDir)

if csv:
    A02.createCsvCorrelationsFromList(name,pdbDir, geos, name, pdbs,250, False)

print('load csv', "Csv/" + name + "_" + name + ".csv")
df_contacts = pd.read_csv("Csv/" + name + "_" + name + ".csv")

if html:
    A03.createHtmlFromPairs(name,name,df_contacts,pairs,geos)





