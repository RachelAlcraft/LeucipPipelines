#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################### run params
CAP = -1
################### set params
ID = 'SG_NZ_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "SG:{NZ}"
GeoB = "SG:{SG,N,NE,NH1,NH2,ND2,NE2,ND1,NZ,NE1,O,OG1,OD1,OD2,OE1,OE2,OG,OH}+1"
GeoAx = "SG:NZ"
GeoBx = "SG:CC"
GeoCx = "NZ:CC"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,4
GeoBMin,GeoBMax = -1,4
GeoCMin,GeoCMax = -1,-1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Overlay = True
HtmlStart = -1
HtmlEnd = -1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# LeucipPlus settings
ExePath ='C:/Dev/Github/PsuMaxima/Linux/out/build/x64-Release/PsuMaxima.exe'
Ccp4Directory = "C:/Dev/Github/ProteinDataFiles/ccp4_data/"
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
LogDirectory = "C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Log/"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import pandas as pd
import A01MakeCsv as A01
import A03RunLeucipPlus as A02
import A04CreateHtml as A03

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
make,run,report = True,True,True

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoAx,GeoBx,GeoCx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")
df_leuci = pd.read_csv("Csv/" + ID + "_04_LeucipPlus.csv")

if run:
    A02.runCreateElectronDensity(df_leuci,ExePath,LogDirectory)

if report:
    A03.runCreateReport(df_leuci, df_geometry,HtmlStart,HtmlEnd,ID,GeoAx,GeoBx,GeoCx,Overlay)





