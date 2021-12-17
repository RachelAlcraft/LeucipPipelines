make,report = True,True
#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
AtomType = 'O'
ID = AtomType + '_ED_CC6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO)+2"
GeoB = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@2)+2"
GeoC = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@3)+2"
GeoD = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@4)+2"
GeoE = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@5)+2"
GeoF = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@6)+2"
GeoAx = AtomType + ":CC1"
GeoBx = AtomType + ":CC2"
GeoCx = AtomType + ":CC3"
GeoDx = AtomType + ":CC4"
GeoEx = AtomType + ":CC5"
GeoFx = AtomType + ":CC6"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,3
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
HB = [4,4,4,4,4,4]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Overlay = True
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

if make:
    A01.runMakeCsv(ID + '',PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoD,GeoE,GeoF,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '_SUMMARY', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)
    A02.runCreateReport(df_geometry,ID+'',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB)


