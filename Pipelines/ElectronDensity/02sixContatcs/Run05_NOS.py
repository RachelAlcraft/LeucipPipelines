make,coords,run,report = True,True,True,True
#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
AtomType = 'SG'
ID = AtomType + '_ED_CC6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_1_1_good.csv'
GeoA = AtomType + ":{NZ}"
GeoB = AtomType + ":(S,O,N,NA,CU,MG,ZN,FE,CO@2)+1"
GeoC = AtomType + ":(S,O,N,NA,CU,MG,ZN,FE,CO@3)+2"
GeoD = AtomType + ":(S,O,N,NA,CU,MG,ZN,FE,CO@4)+3"
GeoE = AtomType + ":(S,O,N,NA,CU,MG,ZN,FE,CO@5)+4"
GeoF = AtomType + ":(S,O,N,NA,CU,MG,ZN,FE,CO@6)+5"
GeoAx = AtomType + ":NZ"
GeoBx = AtomType + ":CC1"
GeoCx = AtomType + ":CC2"
GeoDx = AtomType + ":CC3"
GeoEx = AtomType + ":CC4"
GeoFx = AtomType + ":CC5"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,3
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
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
import A02MakeCoords as A02
import A03RunLeucipPlus as A03
import A04CreateHtml as A04

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoD,GeoE,GeoF,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

if coords:
    A02.runMakeCoords(ID,df_geometry)

df_leuci = pd.read_csv("Csv/" + ID + "_04_LeucipPlus.csv")

if run:
    A03.runCreateElectronDensity(df_leuci,ExePath,LogDirectory)

if report:
    A04.runCreateReport(df_leuci, df_geometry, HtmlStart, HtmlEnd, ID + '', GeoAx,[GeoBx,GeoCx,GeoDx,GeoEx,GeoFx],Overlay)





