#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'NA_GEO6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA2 = ":(FE,S,N,O,NA,CA,CU,MG,ZN)+1"
GeoB2 = ":(FE,S,N,O,NA@2)+1"
GeoC2 = ":(FE,S,N,O,NA@3)+1"
GeoD2 = ":(FE,S,N,O,NA@4)+1"
GeoE2 = ":(FE,S,N,O,NA@5)+1"
GeoF2 = ":(FE,S,N,O,NA@6)+1"
GeoAx = "NA:CC1"
GeoBx = "NA:CC2"
GeoCx = "NA:CC3"
GeoDx = "NA:CC4"
GeoEx = "NA:CC5"
GeoFx = "NA:CC6"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,-1
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
HB = [4,4,4,4,4,4]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
HtmlStart = -1
HtmlEnd = -1
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
make,report = True,True

if make:
    A01.runMakeCsv(ID + '_NA',PdbFile,CAP,PdbDirectory,'(NA)'+GeoA2,'FE'+GeoB2,'FE'+GeoC2,'FE'+GeoD2,'FE'+GeoE2,'FE'+GeoF2,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)
    #A01.runMakeCsv(ID + '_FE1', PdbFile, CAP, PdbDirectory, 'FE1' + GeoA2, 'FE1' + GeoB2, 'FE1' + GeoC2, 'FE1' + GeoD2,'FE1' + GeoE2, 'FE1' + GeoF2, GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, GeoAMin, GeoAMax, GeoBMin,GeoBMax, GeoCMin, GeoCMax, Inclusions, Exclusions)
    #A01.runMakeCsv(ID + '_FE2', PdbFile, CAP, PdbDirectory, 'FE2' + GeoA2, 'FE2' + GeoB2, 'FE2' + GeoC2, 'FE2' + GeoD2,'FE2' + GeoE2, 'FE2' + GeoF2, GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, GeoAMin, GeoAMax, GeoBMin,GeoBMax, GeoCMin, GeoCMax, Inclusions, Exclusions)
    #A01.runMakeCsv(ID + '_FE3', PdbFile, CAP, PdbDirectory, 'FE3' + GeoA2, 'FE3' + GeoB2, 'FE3' + GeoC2, 'FE3' + GeoD2,'FE3' + GeoE2, 'FE3' + GeoF2, GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, GeoAMin, GeoAMax, GeoBMin, GeoBMax, GeoCMin, GeoCMax, Inclusions, Exclusions)
    #A01.runMakeCsv(ID + '_FE4', PdbFile, CAP, PdbDirectory, 'FE4' + GeoA2, 'FE4' + GeoB2, 'FE4' + GeoC2, 'FE4' + GeoD2,'FE4' + GeoE2, 'FE4' + GeoF2, GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, GeoAMin, GeoAMax, GeoBMin,GeoBMax, GeoCMin, GeoCMax, Inclusions, Exclusions)


df_geometry = pd.read_csv("Csv/" + ID + "_NA_03_Geometry.csv")
#df_geometryFE1 = pd.read_csv("Csv/" + ID + "_FE1_03_Geometry.csv")
#df_geometryFE2 = pd.read_csv("Csv/" + ID + "_FE2_03_Geometry.csv")
#df_geometryFE3 = pd.read_csv("Csv/" + ID + "_FE3_03_Geometry.csv")
#df_geometryFE4 = pd.read_csv("Csv/" + ID + "_FE4_03_Geometry.csv")

#df_geometryFE = df_geometryFE[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF']]
#df_geometryFE1 = df_geometryFE1[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF']]
#df_geometryFE2 = df_geometryFE2[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF']]
#df_geometryFE3 = df_geometryFE3[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF']]
#df_geometryFE4 = df_geometryFE4[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF']]

#df_geometry = pd.concat([df_geometryFE,df_geometryFE1,df_geometryFE2,df_geometryFE3,df_geometryFE4])


if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '_SUMMARY', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)
    A02.runCreateReport(df_geometry,ID+'',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB)






