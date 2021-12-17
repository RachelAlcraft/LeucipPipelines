make,coords,run,report = True,True,True,True
#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
AtomType = 'O'
ID = AtomType + '_ED_CC6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@2)"
GeoB = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@3)"
GeoC = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@4)"
GeoD = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@5)"
GeoE = AtomType + ":(S,N,O,NA,CU,MG,ZN,FE,CO@6)"
GeoF = AtomType + ":{OD1,OD2,OE1,OE2}"
GeoAx = AtomType + ":CC1"
GeoBx = AtomType + ":CC2"
GeoCx = AtomType + ":CC3"
GeoDx = AtomType + ":CC4"
GeoEx = AtomType + ":CC5"
GeoFx = AtomType + ":HB"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,3
GeoBMin,GeoBMax = -1,3
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





