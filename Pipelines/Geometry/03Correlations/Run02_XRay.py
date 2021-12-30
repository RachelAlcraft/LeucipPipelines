import pandas as pd

unzip, csv, html  = False, True, True
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
geos = ['N:N+1', 'N:O', 'CB:O', 'N:CA:C', 'N:CA:C:N+1', 'C-1:N:CA:C', 'C-1:C','CA-2:CA-1:CA','CA:CA+1:CA+2']
pairs = [['N:CA:C:N+1','N:N+1'],['N:O','CB:O'],['C-1:N:CA:C','C-1:C'],['CA-2:CA-1:CA','CA:CA+1:CA+2']]
chunk = 200
################################################
pdbDir = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
PdbFile = '../../0ConfigData/Pdbs_1_1_good.csv'
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import A01Unzip as A01
import A02CreateCsv as A02
import A03CreateHtml as A03

dfs = []
pdb_df = pd.read_csv(PdbFile)
pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
pdbs = pdb_df['PDB'].values


if unzip:
    A01.unzipProteome(pdbDir)

if csv:
    A02.createCsvCorrelationsFromList('xray',pdbDir, geos, "Contacts", pdbs,250, False)

print('load csv', "Csv/" + 'xray' + "_Contacts.csv")
df_contacts = pd.read_csv("Csv/" + 'xray' + "_Contacts.csv")
dfs.append(['nmr', df_contacts])

if html:
    A03.createHtmlFromPairs('xray', 'Contacts',df_contacts,pairs,geos)





