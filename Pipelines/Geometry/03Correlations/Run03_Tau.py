import pandas as pd

unzip, csv, html  = False, False, True
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
geoA = 'N:N+1'
geoB = 'C:O'
geoC = 'C:N+1'
geoD = 'N:CA:C'
geoE = 'CA:C:N+1'
geoF = 'N:CA:C:N+1'
geos = [geoA,geoB,geoC,geoD,geoE,geoF]
triples = [[geoA,geoB,geoC],[geoA,geoB,geoD],[geoA,geoB,geoE],[geoA,geoC,geoD],[geoA,geoC,geoE],[geoA,geoD,geoE],[geoB,geoC,geoD],[geoB,geoC,geoE],[geoB,geoD,geoE],[geoC,geoD,geoE],[geoA,geoD,geoF],[geoD,geoC,geoF]]
#                                                                                                                                                                                                              triples = [[geoB,geoA,geoC],[geoB,geoA,geoD],[geoB,geoA,geoE],[geoB,geoA,geoF],[geoB,geoC,geoD],[geoB,geoC,geoE],[geoB,geoD,geoE],[geoB,geoD,geoF],[geoB,geoE,geoF]]
chunk = 500
name = 'tau'
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
    A02.createCsvCorrelationsFromList(name,pdbDir, geos, name, pdbs,250, True)

print('load csv')
df = pd.read_csv("Csv/" + name + "_" + name + ".csv")
df = df.query("bfactor<5")
df = df.query("`C:N+1`<1.5")
df = df.query("`N:N+1`<5")
df = df.query("`C:O`<1.35")



if html:
    A03.createHtmlFromTriples(name, name,df,triples)





