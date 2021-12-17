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
def runMakeCsvMulti(ID,PdbFile,CAP,PdbDirectory,Geos,Inclusions,Exclusions):
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
    geos = Geos
    data = geo.calculateGeometry(geos, log=1)
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)

    # make an id out of the associated atoms
    data['CLASS'] = data.apply(lambda row: help.applyCLASS(pdb_df,row['pdb_code']), axis=1)

    data.to_csv("Csv/" + ID + "_02_Geometry.csv", index=False)

    if Inclusions != {}:
        data = geo.filterDataFrame(data, inclusions=Inclusions)
    if Exclusions != {}:
        data = geo.filterDataFrame(data, exclusions=Exclusions)

    data.to_csv("Csv/" + ID + "_03_Geometry.csv", index=False)

#***************************************************************************************************************
def runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoAx,GeoBx,GeoCx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions):
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
    geos = [GeoA, GeoB,GeoC]
    data = geo.calculateGeometry(geos, log=1)
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)

    data[GeoAx] = data[GeoA]
    data[GeoBx] = data[GeoB]
    data[GeoCx] = data[GeoC]

    # make an id out of the associated atoms
    data['atom'] = data.apply(lambda row: help.applyATOM(row['info' + GeoA], 0), axis=1)
    data['ridA'] = data.apply(lambda row: help.applyRID(row['info' + GeoA]), axis=1)
    data['ridGapA'] = data.apply(lambda row: help.applyRIDGAP(row['info' + GeoA]), axis=1)
    data['aaA'] = data.apply(lambda row: help.applyAA(row['info' + GeoA]), axis=1)
    data['chainA'] = data.apply(lambda row: help.applyCHAIN(row['info' + GeoA]), axis=1)
    data['atomA'] = data.apply(lambda row: help.applyATOM(row['info' + GeoA],1), axis=1)

    data['ridB'] = data.apply(lambda row: help.applyRID(row['info' + GeoB]), axis=1)
    data['ridGapB'] = data.apply(lambda row: help.applyRIDGAP(row['info' + GeoB]), axis=1)
    data['aaB'] = data.apply(lambda row: help.applyAA(row['info' + GeoB]), axis=1)
    data['chainB'] = data.apply(lambda row: help.applyCHAIN(row['info' + GeoB]), axis=1)
    data['atomB'] = data.apply(lambda row: help.applyATOM(row['info' + GeoB],1), axis=1)

    data['ridC'] = data.apply(lambda row: help.applyRID(row['info' + GeoC]), axis=1)
    data['ridGapC'] = data.apply(lambda row: help.applyRIDGAP(row['info' + GeoC]), axis=1)
    data['aaC'] = data.apply(lambda row: help.applyAA(row['info' + GeoC]), axis=1)
    data['chainC'] = data.apply(lambda row: help.applyCHAIN(row['info' + GeoC]), axis=1)
    data['atomC'] = data.apply(lambda row: help.applyATOM(row['info' + GeoC], 1), axis=1)

    data['CLASS'] = data.apply(lambda row: help.applyCLASS(pdb_df,row['pdb_code']), axis=1)

    data['InfoA'] = data['info'+GeoA]
    data['InfoB'] = data['info' + GeoB]
    data['InfoC'] = data['info' + GeoC]

    data.to_csv("Csv/" + ID + "_02_Geometry.csv", index=False)

    if GeoAMin != -1:
        data = data.query('`' + GeoA + '` > ' + str(GeoAMin))
    if GeoAMax != -1:
        data = data.query('`' + GeoA + '` <= ' + str(GeoAMax))
    if GeoBMin != -1:
        data = data.query('`' + GeoB + '` > ' + str(GeoBMin))
    if GeoBMax != -1:
        data = data.query('`' + GeoB + '` <= ' + str(GeoBMax))
    if GeoCMin != -1:
        data = data.query('`' + GeoCx + '` > ' + str(GeoCMin))
    if GeoCMax != -1:
        data = data.query('`' + GeoCx + '` <= ' + str(GeoCMax))
    if Inclusions != {}:
        data = geo.filterDataFrame(data, inclusions=Inclusions)
    if Exclusions != {}:
        data = geo.filterDataFrame(data, exclusions=Exclusions)

    data.to_csv("Csv/" + ID + "_03_Geometry.csv", index=False)

#***************************************************************************************************************
