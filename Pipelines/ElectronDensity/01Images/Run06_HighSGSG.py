#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################### run params
CAP = -1
################### set params
ID = 'SG_SG_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "SG:{SG}+1"
GeoB = "SG:(FE,N,O)+1"
GeoAx = "SG:SG2"
GeoBx = "SG:CC1"
GeoCx = "SG2:CC1"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,3.5
GeoBMin,GeoBMax = -1,-1
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
import A02MakeCoords as A02
import A03RunLeucipPlus as A03
import A04CreateHtml as A04

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
make,coords,run,report = False,True,True,True

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoAx,GeoBx,GeoCx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

if coords:
    A02.runMakeCoords(ID,df_geometry)

df_leuci = pd.read_csv("Csv/" + ID + "_04_LeucipPlus.csv")

if run:
    A03.runCreateElectronDensity(df_leuci,ExePath,LogDirectory)

if report:
    max = len(df_leuci.index)
    chunk = 50
    start = 0
    while start < max:
        xHtmlStart = start
        xHtmlEnd = start + chunk
        if xHtmlEnd > max:
            xHtmlEnd = max
        A04.runCreateReport(df_leuci, df_geometry, xHtmlStart, xHtmlEnd, ID, GeoAx, GeoBx, GeoCx, Overlay)
        start += chunk






