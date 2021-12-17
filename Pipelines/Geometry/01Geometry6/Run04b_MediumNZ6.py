#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'NZ_GEO6_MEDIUM'
PdbFile = '../../0ConfigData/PdbFilesCrossLinkSGNZ_good.csv'
GeoA = "NZ:(FE,S,N,O)+1"
GeoB = "NZ:(FE,S,N,O@2)+1"
GeoC = "NZ:(FE,S,N,O@3)+1"
GeoD = "NZ:(FE,S,N,O@4)+1"
GeoE = "NZ:(FE,S,N,O@5)+1"
GeoF = "NZ:(FE,S,N,O@6)+1"
GeoAx = "NZ:CC1"
GeoBx = "NZ:CC2"
GeoCx = "NZ:CC3"
GeoDx = "NZ:CC4"
GeoEx = "NZ:CC5"
GeoFx = "NZ:CC6"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,3
GeoBMin,GeoBMax = -1,3
GeoCMin,GeoCMax = -1,3
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
    df_geometry = df_geometry.query('`' + GeoA + '` <= ' + str(3.1),hue='aa')
    A02.runCreateReport(df_geometry,ID+'',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB)





