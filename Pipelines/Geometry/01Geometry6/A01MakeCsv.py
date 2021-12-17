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
def runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoD,GeoE,GeoF,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions):
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
    geos = [GeoA, GeoB,GeoC,GeoD,GeoE,GeoF]
    data = geo.calculateGeometry(geos, log=1)
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)

    data[GeoAx] = data[GeoA]
    data[GeoBx] = data[GeoB]
    data[GeoCx] = data[GeoC]
    data[GeoDx] = data[GeoD]
    data[GeoEx] = data[GeoE]
    data[GeoFx] = data[GeoF]

    # make an id out of the associated atoms
    data['atom'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoA], 0), axis=1)

    data['ridA'] = data.apply(lambda row: help.applyDetailFromOther('RID', row['info' + GeoA], 1), axis=1)
    data['ridGapA'] = data.apply(lambda row: help.applyDetailFromOther('GAP', row['info' + GeoA], 1), axis=1)
    data['aaA'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['info' + GeoA], 1), axis=1)
    data['chainA'] = data.apply(lambda row: help.applyDetailFromOther('CHAIN', row['info' + GeoA], 1), axis=1)
    data['atomA'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoA], 1), axis=1)

    data['ridB'] = data.apply(lambda row: help.applyDetailFromOther('RID', row['info' + GeoB], 1), axis=1)
    data['aaB'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['info' + GeoB], 1), axis=1)
    data['chainB'] = data.apply(lambda row: help.applyDetailFromOther('CHAIN', row['info' + GeoB], 1), axis=1)
    data['atomB'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoB], 1), axis=1)

    data['ridC'] = data.apply(lambda row: help.applyDetailFromOther('RID', row['info' + GeoC], 1), axis=1)
    data['aaC'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['info' + GeoC], 1), axis=1)
    data['chainC'] = data.apply(lambda row: help.applyDetailFromOther('CHAIN', row['info' + GeoC], 1), axis=1)
    data['atomC'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoC], 1), axis=1)

    data['ridD'] = data.apply(lambda row: help.applyDetailFromOther('RID', row['info' + GeoD], 1), axis=1)
    data['aaD'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['info' + GeoD], 1), axis=1)
    data['chainD'] = data.apply(lambda row: help.applyDetailFromOther('CHAIN', row['info' + GeoD], 1), axis=1)
    data['atomD'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoD], 1), axis=1)

    data['ridE'] = data.apply(lambda row: help.applyDetailFromOther('RID', row['info' + GeoE], 1), axis=1)
    data['aaE'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['info' + GeoE], 1), axis=1)
    data['chainE'] = data.apply(lambda row: help.applyDetailFromOther('CHAIN', row['info' + GeoE], 1), axis=1)
    data['atomE'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoE], 1), axis=1)

    data['ridF'] = data.apply(lambda row: help.applyDetailFromOther('RID', row['info' + GeoF], 1), axis=1)
    data['aaF'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['info' + GeoF], 1), axis=1)
    data['chainF'] = data.apply(lambda row: help.applyDetailFromOther('CHAIN', row['info' + GeoF], 1), axis=1)
    data['atomF'] = data.apply(lambda row: help.applyDetailFromOther('ATOM', row['info' + GeoF], 1), axis=1)


    data['CLASS'] = data.apply(lambda row: help.applyCLASS(pdb_df,row['pdb_code']), axis=1)

    data['InfoA'] = data['info'+GeoA]
    data['InfoB'] = data['info' + GeoB]
    data['InfoC'] = data['info' + GeoC]
    data['InfoD'] = data['info' + GeoD]
    data['InfoE'] = data['info' + GeoE]
    data['InfoF'] = data['info' + GeoF]

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
