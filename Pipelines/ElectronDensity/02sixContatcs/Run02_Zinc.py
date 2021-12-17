#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'ZN_ED_CC6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "(ZN):(S,N,O,NA,CU,MG,ZN,FE@2)"
GeoB = "(ZN):(S,N,O,NA,CU,MG,ZN,FE@3)"
GeoC = "(ZN):(S,N,O,NA,CU,MG,ZN,FE@4)"
GeoD = "(ZN):(S,N,O,NA,CU,MG,ZN,FE@5)"
GeoE = "(ZN):(S,N,O,NA,CU,MG,ZN,FE@6)"
GeoF = "(ZN):(S,N,O,NA,CU,MG,ZN,FE@7)"
GeoAx = "ZN:CC1"
GeoBx = "ZN:CC2"
GeoCx = "ZN:CC3"
GeoDx = "ZN:CC4"
GeoEx = "ZN:CC5"
GeoFx = "ZN:CC6"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,3
GeoBMin,GeoBMax = -1,3
GeoCMin,GeoCMax = -1,-1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Overlay = True
HtmlStart = -1
HtmlEnd = 100
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
make,coords,run,report = False,False,False,True

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





