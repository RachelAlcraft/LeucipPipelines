'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import sys
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Helpers as help
import Bio.PDB as bio
import pandas as pd
from LeucipPy import GeoDataFrame as gdf
from LeucipPy import LeucipPy as leu

#***************************************************************************************************************
def runMakeCsv(ID,PdbFile,PdbDir,AllGeos,aa_filter=[]):
    print('LeucipPipeline: - Creating data -----------------------------------------------------------')
    pdb_df = pd.read_csv(PdbFile)
    pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
    pdbs = pdb_df['PDB'].values
    #pdbs = pdbs[:10]
    ### 2) Loop through and create bio python objects
    parser = bio.PDBParser()
    strucs = []
    count = 1
    for pdb in pdbs:
        pdb = pdb.lower()
        print("CL: creating structures", count, '/', len(pdbs))
        count += 1
        pdb_file, pdb_html_loc = leu.getPdbLink(pdb)
        struc = parser.get_structure(pdb, PdbDir + pdb_file)
        strucs.append(struc)

    geo = gdf.GeoDataFrame(strucs, log=0)
    data = geo.calculateGeometry(AllGeos, log=1)

    if len(aa_filter) > 0:
        print('filter on',aa_filter)
        data = data[data['aa'].isin(aa_filter)]

    print('Save to', "Csv/" + ID + "_01_Geometry.csv")
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)


#***************************************************************************************************************
def runMakeGeosPairwise(atoms_list):
    #Make a pairwise listout of atoms
    geos = []
    for i in range(0, len(atoms_list)):
        atmA = atoms_list[i]
        for j in range(i+1, len(atoms_list)):
            atmB = atoms_list[j]
            if atmA != atmB:
                geos.append(atmA + ':' + atmB)
    return geos
