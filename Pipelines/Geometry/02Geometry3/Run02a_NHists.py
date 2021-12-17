#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ID = 'N_HISTS_HIGH'
PdbFile = '../../0ConfigData/Pdbs_Under1_good.csv'
Geos = ['N:O-7','N:O-6','N:O-5','N:O-4','N:O-3','N:O-2','N:O-1','N:O','N:O+1','N:O+2','N:O+3','N:O+4','N:O+5','N:O+6','N:O+7']

Inclusions = {}
Exclusions = {}
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
    A01.runMakeCsvMulti(ID,PdbFile,CAP,PdbDirectory,Geos,Inclusions,Exclusions)

df_geometry = pd.read_csv("Csv/" + ID + "_03_Geometry.csv")

# We need to remove the contacts so far out that they don;t exists, ie chain ends

#for i in range(0,len(Geos)):
#    df_geometry = df_geometry.query('`' + Geos[i] + '` < 10')
#df_geometry = df_geometry.query('`N:O-4` < 6')
#df_geometry = df_geometry.query('`N:O-3` < 6')
#df_geometry = df_geometry.query('`N:O-1` < 3')
#df_geometry = df_geometry.query('`N:O` < 5')
#df_geometry = df_geometry.query('`N:O+1` < 10')


if report:
    A02.runCreateReportAnalysisMulti(df_geometry, ID + '_MULTI', Geos)






