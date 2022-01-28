'''
'''
import shutil

PdbFile = '../../0ConfigData/Pdbs_1_1_good.csv'
OutDir = 'C:/Dev/Github/ProteinDataFiles/pdb_alpha/HighRes/'
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
####################################################################################
from LeucipPy import LeucipPy as leu
import pandas as pd
import os
from urllib.request import urlretrieve

pdbdic = {}

pdb_df = pd.read_csv(PdbFile)
pdbs =pdb_df['PDB'].values
prots =pdb_df['CLASS'].values

### 3) Loop through and ensure files are retrieved

good_pdbs = []
for count in range(0,len(pdbs)):
#for count in range(0, 10):
    pdb = pdbs[count]
    print("LeucipPipelines:Utilities:01Download(1) Downloading",pdb,count, '/', len(pdbs))
    prot = prots[count]
    pdb = pdb.lower()
    pdb_file, pdb_html_loc = leu.getPdbLink(pdb)
    #print(pdb_file,pdb_html_loc)
    pdb_exists = False

    if not os.path.exists(PdbDirectory + pdb_file):
        try:
            urlretrieve(pdb_html_loc, PdbDirectory + pdb_file)
            pdb_exists = True
        except:
            print("...!!! No PDB data for", pdb)
    else:
        pdb_exists = True

    if pdb_exists:
        print("...adding", pdb)
        shutil.copy(PdbDirectory + pdb_file, OutDir + pdb_file)



