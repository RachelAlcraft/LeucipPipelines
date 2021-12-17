#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'SD_GEO6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "SD:(FE,S,N,O)+1"
GeoB = "SD:(FE,S,N,O@2)+1"
GeoC = "SD:(FE,S,N,O@3)+1"
GeoD = "SD:(FE,S,N,O@4)+1"
GeoE = "SD:(FE,S,N,O@5)+1"
GeoF = "SD:(FE,S,N,O@6)+1"
GeoAx = "SD:CC1"
GeoBx = "SD:CC2"
GeoCx = "SD:CC3"
GeoDx = "SD:CC4"
GeoEx = "SD:CC5"
GeoFx = "SD:CC6"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,-1
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
HB = [3,3,3,3,3,3]
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
make,report = False,True

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoD,GeoE,GeoF,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

if report:
    import pandas as pd

    A02.runCreateReportAnalysis(df_geometry, ID + '_SUMMARY', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)
    df_geometry = df_geometry.query('`' + GeoAx + '` <= ' + str(3.1))
    A02.runCreateReport(df_geometry,ID+'',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB)





