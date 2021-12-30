#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'GLN_GEO6_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
GeoA2 = ":(S,N,O)+2"
GeoB2 = ":(S,N,O@2)+2"
GeoC2 = ":(S,N,O@3)+2"
GeoD2 = ":(S,N,O@4)+2"
GeoE2 = ":(S,N,O@5)+2"
GeoF2 = ":(S,N,O@6)+2"
GeoAx = "ON:CC1"
GeoBx = "ON:CC2"
GeoCx = "ON:CC3"
GeoDx = "ON:CC4"
GeoEx = "ON:CC5"
GeoFx = "ON:CC6"
Inclusions = {'aa':['GLN']}
Exclusions = {}
GeoAMin, GeoAMax = -1,-1
GeoBMin,GeoBMax = -1,-1
GeoCMin,GeoCMax = -1,-1
HB = [3,3,3,3,3,3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
make,report = False,True

if make:
    A01.runMakeCsv(ID + '_OE1',PdbFile,CAP,PdbDirectory,'OE1'+GeoA2,'OE1'+GeoB2,'OE1'+GeoC2,'OE1'+GeoD2,'OE1'+GeoE2,'OE1'+GeoF2,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)
    A01.runMakeCsv(ID + '_NE2',PdbFile,CAP,PdbDirectory,'NE2'+GeoA2,'NE2'+GeoB2,'NE2'+GeoC2,'NE2'+GeoD2,'NE2'+GeoE2,'NE2'+GeoF2,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,GeoAMin,GeoAMax,GeoBMin,GeoBMax,GeoCMin,GeoCMax,Inclusions,Exclusions)

df_geometryOE1 = pd.read_csv("Csv/" + ID + "_OE1_03_Geometry.csv")
df_geometryNE2 = pd.read_csv("Csv/" + ID + "_NE2_03_Geometry.csv")

df_geometryOE1 = df_geometryOE1[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF','ridGapA']]
df_geometryNE2 = df_geometryNE2[['pdb_code','resolution','chain','aa','rid','bfactor','occupancy',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,'atom','ridA','aaA','chainA','atomA','ridB','aaB','chainB','atomB','ridC','aaC','chainC','atomC','ridD','aaD','chainD','atomD','ridE','aaE','chainE','atomE','ridF','aaF','chainF','atomF','CLASS','InfoA','InfoB','InfoC','InfoD','InfoE','InfoF','ridGapA']]

df_geometry = pd.concat([df_geometryOE1,df_geometryNE2])
df_geometry.to_csv("Csv/" + ID + "_03_Geometry.csv", index=False)

# now we want to make 2 reports, one where there are 2 close contacts and one where there are not
df_geometry2 = df_geometry.query('`' + GeoAx + '` < ' + str(3))
#df_geometry2 = df_geometry2.query('`' + GeoBx + '` < ' + str(3))

#df_geometry0 = df_geometry.query('`' + GeoBx + '` >= ' + str(3))

if report:
    A02.runCreateReportAnalysis(df_geometry, ID + '', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)
    A02.runCreateReport(df_geometry2,ID+'_2C',GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB)
    A02.runCreateReport(df_geometry, ID + '_All', GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB)







