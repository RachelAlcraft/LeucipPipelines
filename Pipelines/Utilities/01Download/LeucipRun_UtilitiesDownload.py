'''


'''

import Config as cfg
from LeucipPy import LeucipPy as leu
import pandas as pd
import os
from urllib.request import urlretrieve

pdbdic = {}

pdb_df = pd.read_csv(cfg.PdbFile)
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
    ccp4_exists = False
    diff_exists = False
    if not os.path.exists(cfg.PdbDirectory + pdb_file):
        try:
            urlretrieve(pdb_html_loc, cfg.PdbDirectory + pdb_file)
            pdb_exists = True
        except:
            print("...!!! No PDB data for", pdb)
    else:
        pdb_exists = True


    ccp4_file, ccp4_html_loc = leu.getElectronDensityLink(pdb)
    if not os.path.exists(cfg.Ccp4Directory + ccp4_file):
        try:
            urlretrieve(ccp4_html_loc, cfg.Ccp4Directory + ccp4_file)
            ccp4_exists = True
        except:
            print("...No ccp4 data for", pdb)
    else:
        ccp4_exists = True

    ccp4_difffile, ccp4_diffhtml_loc = leu.getElectronDensityLink(pdb,diff=True)
    if not os.path.exists(cfg.Ccp4Directory + ccp4_difffile):
        try:
            urlretrieve(ccp4_diffhtml_loc, cfg.Ccp4Directory + ccp4_difffile)
            diff_exists = True
        except:
            print("...No ccp4 diff data for", pdb)
    else:
        diff_exists = True

    if pdb_exists and ccp4_exists and diff_exists:
        print("...adding", pdb)
        good_pdbs.append([pdb,prot])

df_good = pd.DataFrame(good_pdbs, columns=['PDB','CLASS'])
df_good.to_csv(cfg.OutFile,index=False)
print("Saved to",cfg.OutFile)
