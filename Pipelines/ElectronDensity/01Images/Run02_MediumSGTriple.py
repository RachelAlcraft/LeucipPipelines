#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################### run params
CAP = -1
################### set params
ID = 'SG_TRIPLE_MEDIUM'
PdbFile = '../../0ConfigData/PdbFilesCrossLinkSGNZ_good.csv'
GeoA = "SG:{SG,NZ,N,NE,NH1,NH2,ND2,NE2,ND1,NZ,NE1,O,OG1,OD1,OD2,OE1,OE2,OG,OH}+1"
GeoB = "SG:{SG,NZ,N,NE,NH1,NH2,ND2,NE2,ND1,NZ,NE1,O,OG1,OD1,OD2,OE1,OE2,OG,OH@2}+1"
GeoAx = "SG:CC1"
GeoBx = "SG:CC2"
GeoCx = "CC1:CC2"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,4
GeoBMin,GeoBMax = -1,4
GeoCMin,GeoCMax = -1,3.2
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

make,run,report = False,False,True

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoAx,GeoBx,GeoCx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")
df_leuci = pd.read_csv("Csv/" + ID + "_04_LeucipPlus.csv")

if run:
    A02.runCreateElectronDensity(df_leuci,ExePath,LogDirectory)

if report:
    max = len(df_leuci.index)
    chunk = 50
    start = 250
    while start < max:
        xHtmlStart = start
        xHtmlEnd = start+chunk
        if xHtmlEnd > max:
            xHtmlEnd = max
        A03.runCreateReport(df_leuci, df_geometry,xHtmlStart,xHtmlEnd,ID,GeoAx,GeoBx,GeoCx,Overlay)
        start += chunk



