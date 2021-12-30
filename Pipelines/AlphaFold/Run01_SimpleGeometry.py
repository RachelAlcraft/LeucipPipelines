import gc
import sys
import pandas as pd
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Log as log

unzip,csv,html,compare = False,False,True,False
csv_rewrite = False
chunk = 500
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#names = ['EColi','Yeast','Drome','Mouse','Human']
names = ['pdb_data','nmr_data']
geosA = ['C:N+1', 'C:O', 'CA:C', 'N:CA', 'N:CA:C', 'N:CA:C:N+1', 'C-1:N:CA:C', 'CA:C:N+1']
geosB = ['N:N+1', 'N:O', 'CB:O', 'N:CA:C', 'N:CA:C:N+1', 'C-1:N:CA:C', 'C-1:C','CA-2:CA-1:CA','CA:CA+1:CA+2']
geoChecksA = [['C:N+1',1,2], ['C:O',1,2], ['CA:C',1,2], ['N:CA',1,2], ['N:CA:C',90,130]]
geoChecksB = [['CA-2:CA-1:CA',60,170],['CA:CA+1:CA+2',60,170],['N:O',2,4], ['CB:O',2,4]]
################################################
#dir = 'C:/Dev/Github/ProteinDataFiles/pdb_alpha/'
dir = 'C:/Dev/Github/ProteinDataFiles/'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import A01Unzip as A01
import A02CreateCsv as A02
import A03CreateHtml as A03

dfs = []
for name in names:
    print(log.getTime(),name)

    if unzip:
        A01.unzipProteome(name,dir)

    if csv:
        A02.createCsvCorrelations(name,dir,geosA,"Geometry",chunk,csv_rewrite)
        gc.collect(generation=2)
        A02.createCsvCorrelations(name, dir, geosB,"Correlations",chunk, csv_rewrite)
        gc.collect(generation=2)

    print(log.getTime(),'load csv',"Csv/" + name + "_Geometry.csv")
    df_geometry = pd.read_csv("Csv/" + name + "_Geometry.csv")
    print(log.getTime(),'load csv', "Csv/" + name + "_Correlations.csv")
    df_correlation = pd.read_csv("Csv/" + name + "_Correlations.csv")
    dfs.append([name, df_geometry,df_correlation])


    if html:
        print(log.getTime(),'create html',name)
        df_geometry['Probability'] = 100-df_geometry['bfactor']
        df_correlation['Probability'] = 100-df_correlation['bfactor']
        A03.createHtml(name,df_geometry,df_correlation,geoChecksA,geoChecksB)

if compare:
    print(log.getTime(),'Start HTML compare')
    A03.compareProteomes(dfs)


    
    
