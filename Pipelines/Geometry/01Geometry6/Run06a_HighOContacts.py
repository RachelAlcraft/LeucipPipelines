#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'O_WITH_GEO6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "O:(FE,S,N,O)+2"
GeoB = "O:(FE,S,N,O@2)+2"
GeoC = "O:(FE,S,N,O@3)+2"
GeoD = "O:(FE,S,N,O@4)+2"
GeoE = "O:(FE,S,N,O@5)+2"
GeoF = "O:(FE,S,N,O@6)+2"
GeoAx = "O:CC1"
GeoBx = "O:CC2"
GeoCx = "O:CC3"
GeoDx = "O:CC4"
GeoEx = "O:CC5"
GeoFx = "O:CC6"
Inclusions = {}
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
make,report = True,True

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoD,GeoE,GeoF,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

# now we want to make 2 reports, one where there are 2 close contacts and one where there are not
df_geometry = df_geometry.query("aa != 'HOH'")
df_geometry2 = df_geometry.query('`' + GeoAx + '` < ' + str(3))
df_geometry2 = df_geometry2.query('`' + GeoBx + '` < ' + str(3))

if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '_SUMMARY', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)
    A02.runCreateReport(df_geometry2,ID+'_2C',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB)





