#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'NO_TRIPLE_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "N:O+1"
GeoB = "N:O-1"
GeoC = "N:(O)+2"
GeoAx = "N:OP"
GeoBx = "N:OM"
GeoCx = "N:OC"
Inclusions = {}
Exclusions = {}
GeoAMin, GeoAMax = -1,-1
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
HB = [3,3,3,3,3,3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
CAP = -1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import pandas as pd
import A01MakeCsv as A01
import A02CreateHtml as A02
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
make,report = False,True

if make:
    A01.runMakeCsv(ID,PdbFile,CAP,PdbDirectory,GeoA,GeoB,GeoC,GeoAx,GeoBx,GeoCx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

# We need to remove the contacts so far out that they don;t exists, ie chain ends
df_geometry = df_geometry.query('`' + GeoAx + '` < ' + str(8))
df_geometry = df_geometry.query('`' + GeoBx + '` < ' + str(6))

if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '_SUMMARY', GeoAx, GeoBx, GeoCx, HB)
    A02.runCreateReport(df_geometry, ID + '', GeoAx, GeoBx, GeoCx, HB)






