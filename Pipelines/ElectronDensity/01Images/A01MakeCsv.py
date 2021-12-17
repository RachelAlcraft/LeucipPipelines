'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import Helpers as help
import Bio.PDB as bio
import pandas as pd
from LeucipPy import GeoDataFrame as gdf
from LeucipPy import LeucipPy as leu

#***************************************************************************************************************
def runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoAx,GeoBx,GeoCx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions):
    print('LeucipPipeline:01Images - 1) Creating data -----------------------------------------------------------')
    pdb_df = pd.read_csv(PdbFile)
    pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
    pdbs = pdb_df['PDB'].values
    ### 2) Loop through and create bio python objects
    parser = bio.PDBParser()
    strucs = []
    count = 1
    if CAP!=-1:#while testing reduce data
        pdbs = pdbs[:CAP]
    for pdb in pdbs:
        pdb = pdb.lower()
        print("CL: creating structures", count, '/', len(pdbs))
        count += 1
        pdb_file, pdb_html_loc = leu.getPdbLink(pdb)
        struc = parser.get_structure(pdb, PdbDirectory + pdb_file)
        strucs.append(struc)

    geo = gdf.GeoDataFrame(strucs, log=0)
    geos = [GeoA, GeoB]
    data = geo.calculateGeometry(geos, log=1)
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)

    data[GeoAx] = data[GeoA]
    data[GeoBx] = data[GeoB]

    # make an id out of the associated atoms
    data['atom1'] = data.apply(lambda row: help.applyATOM(row['info' + GeoA], 0), axis=1)
    data['rid2'] = data.apply(lambda row: help.applyRID(row['info' + GeoA]), axis=1)
    data['aa2'] = data.apply(lambda row: help.applyAA(row['info' + GeoA]), axis=1)
    data['chain2'] = data.apply(lambda row: help.applyCHAIN(row['info' + GeoA]), axis=1)
    data['atom2'] = data.apply(lambda row: help.applyATOM(row['info' + GeoA],1), axis=1)

    data['rid3'] = data.apply(lambda row: help.applyRID(row['info' + GeoB]), axis=1)
    data['aa3'] = data.apply(lambda row: help.applyAA(row['info' + GeoB]), axis=1)
    data['chain3'] = data.apply(lambda row: help.applyCHAIN(row['info' + GeoB]), axis=1)
    data['atom3'] = data.apply(lambda row: help.applyATOM(row['info' + GeoB],1), axis=1)

    #data['ridDiff'] = data.apply(lambda row: help.applyRIDDIFF(row['rid2'],row['rid3']), axis=1)
    data['CLASS'] = data.apply(lambda row: help.applyCLASS(pdb_df,row['pdb_code']), axis=1)
    data['ID'] = data.apply(lambda row: help.applyID(row['aa'],row['rid'],row['chain'],row['atom1'],row['aa2'],row['rid2'],row['chain2'],row['atom2'],row['aa3'],row['rid3'],row['chain3'],row['atom3']), axis=1)

    data.to_csv("Csv/" + ID + "_02_Geometry.csv", index=False)

    if GeoAMin != -1:
        data = data.query('`' + GeoA + '` >= ' + str(GeoAMin))
    if GeoAMax != -1:
        data = data.query('`' + GeoA + '` <= ' + str(GeoAMax))
    if GeoBMin != -1:
        data = data.query('`' + GeoB + '` >= ' + str(GeoBMin))
    if GeoBMax != -1:
        data = data.query('`' + GeoB + '` <= ' + str(GeoBMax))
    if Inclusions != {}:
        data = geo.filterDataFrame(data, inclusions=Inclusions)
    if Exclusions != {}:
        data = geo.filterDataFrame(data, exclusions=Exclusions)

    atomsdata = geo.calculateData(log=1)
    data['Coords1'] = data.apply(lambda row: help.applyCOORDS(geo, atomsdata, row['pdb_code'], row['info' + GeoA], 0),axis=1)
    data['Coords2'] = data.apply(lambda row: help.applyCOORDS(geo, atomsdata, row['pdb_code'], row['info' + GeoA], 1),axis=1)
    data['Coords3'] = data.apply(lambda row: help.applyCOORDS(geo, atomsdata, row['pdb_code'], row['info' + GeoB], 1),axis=1)
    data[GeoCx] = data.apply(lambda row: help.applyDISTANCE(row['Coords2'], row['Coords3']), axis=1)

    if GeoCMin != -1:
        data = data.query('`' + GeoCx + '` >= ' + str(GeoCMin))
    if GeoCMax != -1:
        data = data.query('`' + GeoCx + '` <= ' + str(GeoCMax))

    data.to_csv("Csv/" + ID + "_03_Geometry.csv", index=False)
    atomsdata.to_csv("Csv/" + ID + "_03_Atoms.csv", index=False)

#***************************************************************************************************************
