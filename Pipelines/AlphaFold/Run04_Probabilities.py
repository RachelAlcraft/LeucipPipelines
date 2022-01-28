import gc
import sys
import pandas as pd

sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Log as log

prob = True
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#names = ['EColi','Yeast','Drome','Mouse','Human']
names = ['Human']
#names = ['HighRes']
################################################
scattersA = [['C-1:N:CA:C','N:CA:C:N+1'],['N:CA:C','N:CA:C:N+1']]
scattersB = [['N:O','CB:O'],['CA-2:CA-1:CA','CA:CA+1:CA+2']]
################################################
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_alpha/'
#dir = 'C:/Dev/Github/ProteinDataFiles/'
#dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data_redo/'
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import A03CreateHtml as A03


for name in names:
    print(log.getTime(), name)

    print(log.getTime(), 'load csv', "Csv/" + name + "_Geometry.csv")
    df_geometry = pd.read_csv("Csv/" + name + "_Geometry.csv")
    df_geometryA = df_geometry.query("aa=='GLY'")
    #df_geometryB = df_geometry.query("aa!='GLY'")

    print(log.getTime(), 'load csv', "Csv/" + name + "_Correlations.csv")
    df_correlation = pd.read_csv("Csv/" + name + "_Correlations.csv")
    df_correlationA = df_correlation.query("aa=='GLY'")
    #df_correlationB = df_correlation.query("aa!='GLY'")

    if prob:
        print(log.getTime(), 'create html', name)
        #A03.createProbability(name, 'prob_nogly', df_geometryB, df_correlation, scattersA, scattersB,allone=True)
        A03.createProbability(name, 'prob_gly',df_geometryA, df_correlation,scattersA, scattersB,allone=True)
