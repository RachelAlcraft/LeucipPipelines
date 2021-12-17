#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################### run params
CAP = -1
################### set params
ID = 'SG_HB_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "SG:{SG,NZ,N,NE,NH1,NH2,ND2,NE2,ND1,NZ,NE1,O,OG1,OD1,OD2,OE1,OE2,OG,OH}+1"
GeoB = "SG:{SG,NZ,N,NE,NH1,NH2,ND2,NE2,ND1,NZ,NE1,O,OG1,OD1,OD2,OE1,OE2,OG,OH@2}+1"
GeoAx = "SG:CC1"
GeoBx = "SG:CC2"
GeoCx = "CC1:CC2"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = 3.1,3.8
GeoBMin,GeoBMax = -1,4
GeoCMin,GeoCMax = -1,3
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Overlay = True
HtmlStart = -1
HtmlEnd = 5
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
    A04.runCreateReport(df_leuci, df_geometry,HtmlStart,HtmlEnd,ID+'',GeoAx,GeoBx,GeoCx,Overlay)

'''
if report:
    df_geometryA = df_geometry.query("ridDiff == 1")
    df_geometryB = df_geometry.query("ridDiff == -1")
    df_geometryC = df_geometry.query("ridDiff == 0")
    df_geometryD = df_geometry.query("ridDiff != 1")
    df_geometryD = df_geometryD.query("ridDiff != -1")
    df_geometryD = df_geometryD.query("ridDiff != 0")

    A03.runCreateReport(df_leuci, df_geometryD, HtmlStart, HtmlEnd, ID + '_Other', GeoAx, GeoBx, GeoCx, Overlay)
    A03.runCreateReport(df_leuci, df_geometryA,HtmlStart,HtmlEnd,ID+'_RID1',GeoAx,GeoBx,GeoCx,Overlay)
    A03.runCreateReport(df_leuci, df_geometryB, HtmlStart, HtmlEnd, ID+'_RID_1', GeoAx, GeoBx, GeoCx, Overlay)
    A03.runCreateReport(df_leuci, df_geometryC, HtmlStart, HtmlEnd, ID + '_RID0', GeoAx, GeoBx, GeoCx, Overlay)

'''



