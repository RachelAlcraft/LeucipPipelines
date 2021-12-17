#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'ON_TRIPLE_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA = "O:N+1"
GeoB = "O:N-1"
GeoC = "O:(N)+2"
GeoAx = "O:NP"
GeoBx = "O:NM"
GeoCx = "O:NC"
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

# now we want to make 2 reports, one where there are 2 close contacts and one where there are not
df_geometry = df_geometry.query("aa != 'HOH'")
df_geometry2 = df_geometry.query('`' + GeoAx + '` < ' + str(HB[0]))
df_geometry2 = df_geometry2.query('`' + GeoBx + '` < ' + str(HB[1]))

if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '_SUMMARY', GeoAx, GeoBx, GeoCx, HB)
    A02.runCreateReport(df_geometry, ID + '', GeoAx, GeoBx, GeoCx, HB)
    A02.runCreateReport(df_geometry2,ID+'_2C',GeoAx,GeoBx,GeoCx,HB)





