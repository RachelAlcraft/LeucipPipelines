import gc

import pandas as pd

unzip, csv, html, compare = False, False, False, True
csv_rewrite = False
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
names = ['EColi', 'Yeast', 'Drome', 'Mouse', 'Human']
#names = ['Problems','EColi']
geosA = ['C:N+1', 'C:O', 'CA:C', 'N:CA', 'N:CA:C', 'N:CA:C:N+1', 'C-1:N:CA:C', 'CA:C:N+1']
geosB = ['N:N+1', 'N:O', 'CB:O', 'N:CA:C', 'N:CA:C:N+1', 'C-1:N:CA:C', 'C-1:C', 'CA-2:CA-1:CA', 'CA:CA+1:CA+2']
geoChecksA = [['C:N+1', 1, 2], ['C:O', 1, 2], ['CA:C', 1, 2], ['N:CA', 1, 2], ['N:CA:C', 90, 130]]
geoChecksB = [['CA-2:CA-1:CA', 60, 170], ['CA:CA+1:CA+2', 60, 170], ['N:O', 2, 4], ['CB:O', 2, 4]]
################################################
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_alpha/'
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import A01Unzip as A01
import A02CreateCsv as A02
import A03CreateHtml as A03

dfs = []
for name in names:

    if unzip:
        A01.unzipProteome(name, dir)

    if csv:
        A02.createCsvCorrelations(name, dir, geosA, "Geometry", 500, csv_rewrite)
        gc.collect(generation=2)
        A02.createCsvCorrelations(name, dir, geosB, "Correlations", 500, csv_rewrite)
        gc.collect(generation=2)

    print('load csv', "Csv/" + name + "_Geometry.csv")
    df_geometry = pd.read_csv("Csv/" + name + "_Geometry.csv")
    print('load csv', "Csv/" + name + "_Correlations.csv")
    df_correlation = pd.read_csv("Csv/" + name + "_Correlations.csv")
    dfs.append([name, df_geometry, df_correlation])

    if html:
        A03.createHtml(name, df_geometry, df_correlation, geoChecksA, geoChecksB)

if compare:
    print('Start HTML compare')
    A03.compareProteomes(dfs)




