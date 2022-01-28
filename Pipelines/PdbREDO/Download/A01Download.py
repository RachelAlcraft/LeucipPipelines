import os
import ssl
import pandas as pd
import urllib.request as urlrq
import certifi

def download(pdbDir,PdbFile, OutFile):

    '''
    wget https://pdb-redo.eu/db/9xyz/9xyz_final.pdb downloads the fully optimised (re-refined and rebuilt) PDB file.
    wget https://pdb-redo.eu/db/9xyz/9xyz_final.mtz downloads the MTZ files to generate electron density maps.
    https://pdb-redo.eu/download-info.html;jsessionid=9226A449595FF07E2DCA587EA67F5A16
    '''

    pdb_df = pd.read_csv(PdbFile)
    pdbs = pdb_df['PDB'].values
    prots = pdb_df['CLASS'].values
    redos = []

    for count in range(0, len(pdbs)):
        # for count in range(0, 10):
        pdb = pdbs[count]
        prot = prots[count]
        pdb = pdb.lower()
        print("LeucipPipelines:Utilities:01Download(1) Downloading", pdb, count, '/', len(pdbs))
        pdb_fullfile = pdbDir + "pdb" + pdb + ".ent"
        web_fullfile = "https://pdb-redo.eu/db/" + pdb + "/" + pdb + "_final.pdb"
        print(web_fullfile,pdb_fullfile)
        if not os.path.exists(pdb_fullfile):
            try:
                resp = urlrq.urlopen(web_fullfile,cafile=certifi.where())
                with open(pdb_fullfile,'wb') as output:
                    output.write(resp.read())
                pdb_exists = True
            except:
                print("...!!! No PDB data for", pdb)
                pdb_exists = False
        else:
            pdb_exists = True

        mtz_fullfile = pdbDir + "pdb" + pdb + ".mtz"
        webmtz_fullfile = "https://pdb-redo.eu/db/" + pdb + "/" + pdb + "_final.mtz"
        if not os.path.exists(mtz_fullfile):
            try:
                resp = urlrq.urlopen(webmtz_fullfile, cafile=certifi.where())
                with open(mtz_fullfile, 'wb') as output:
                    output.write(resp.read())
                pdb_exists = pdb_exists and True

            except:
                pdb_exists = False
                print("...!!! No PDB data for", pdb)
        else:
            pdb_exists = pdb_exists and True



        if pdb_exists:
            redos.append([pdb,prot])

    print(redos)
    df_good = pd.DataFrame(redos, columns=['PDB', 'CLASS'])
    df_good.to_csv(OutFile, index=False)
    print("Saved to", OutFile)

