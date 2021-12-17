make,report = True,True
#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
aa = 'ASP'
atoms = ['OD1','OD2']
####################################################################
ID = aa + '_GEO6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA2 = ":(S,N,O,NA,CU,MG,ZN,FE,CO,CA)+2"
GeoB2 = ":(S,N,O,NA,CU,MG,ZN,FE,CO,CA@2)+2"
GeoC2 = ":(S,N,O,NA,CU,MG,ZN,FE,CO,CA@3)+2"
GeoD2 = ":(S,N,O,NA,CU,MG,ZN,FE,CO,CA@4)+2"
GeoE2 = ":(S,N,O,NA,CU,MG,ZN,FE,CO,CA@5)+2"
GeoF2 = ":(S,N,O,NA,CU,MG,ZN,FE,CO,CA@6)+2"
GeoAx = "OD:CC1"
GeoBx = "OD:CC2"
GeoCx = "OD:CC3"
GeoDx = "OD:CC4"
GeoEx = "OD:CC5"
GeoFx = "OD:CC6"
Inclusions = {'aa':[aa]}
Exclusions = {}
GeoAMin, GeoAMax = -1,-1
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
HB = [3,3,3,3,3,3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
CAP = -1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ExePath ='C:/Dev/Github/PsuMaxima/Linux/out/build/x64-Release/PsuMaxima.exe'
Ccp4Directory = "C:/Dev/Github/ProteinDataFiles/ccp4_data/"
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
LogDirectory = "C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Log/"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import pandas as pd
import A01MakeCsv as A01
import A02CreateHtml as A02

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if make:
    for atom in atoms:
        A01.runMakeCsv(ID + '_' + atom, PdbFile, CAP, PdbDirectory, atom + GeoA2, atom + GeoB2, atom + GeoC2,
                       atom + GeoD2, atom + GeoE2, atom + GeoF2, GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, GeoAMin,
                       GeoAMax, GeoBMin, GeoBMax, GeoCMin, GeoCMax, Inclusions, Exclusions)

    dfs = []
    if len(atoms) == 1:
        df_geometryAt = pd.read_csv("Csv/" + ID + "_" + atoms[0] + "_03_Geometry.csv")
        df_geometryAt.to_csv("Csv/" + ID + "_03_Geometry.csv", index=False)

    else:
        for atom in atoms:
            df_geometryAt = pd.read_csv("Csv/" + ID + "_" + atom + "_03_Geometry.csv")
            df_geometryAt = df_geometryAt[
                ['pdb_code', 'resolution', 'chain', 'aa', 'rid', 'bfactor', 'occupancy', GeoAx, GeoBx, GeoCx, GeoDx,
                 GeoEx, GeoFx, 'atom', 'ridA', 'aaA', 'chainA', 'atomA', 'ridB', 'aaB', 'chainB', 'atomB', 'ridC',
                 'aaC', 'chainC', 'atomC', 'ridD', 'aaD', 'chainD', 'atomD', 'ridE', 'aaE', 'chainE', 'atomE', 'ridF',
                 'aaF', 'chainF', 'atomF', 'CLASS', 'InfoA', 'InfoB', 'InfoC', 'InfoD', 'InfoE', 'InfoF', 'ridGapA']]
            dfs.append(df_geometryAt)

        df_geometry = pd.concat(dfs)
        df_geometry.to_csv("Csv/" + ID + "_03_Geometry.csv", index=False)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")
df_geometry2 = df_geometry.query('`' + GeoAx + '` < ' + str(HB[0]))
df_geometry2 = df_geometry2.query('`' + GeoBx + '` < ' + str(HB[1]))

if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)
    A02.runCreateReport(df_geometry2, ID + '_2C', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)