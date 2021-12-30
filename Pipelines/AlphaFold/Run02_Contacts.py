import gc
import sys
import pandas as pd
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Log as log


unzip, csv, html  = False, True, True
overwrite=True
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
names = ['EColi','Yeast','Drome']
geoA = 'N:{N}+2'
geoB = 'O:{O}+2'
geoC = 'N:(O)+2'
chunk = 500
################################################
geos = [geoA,geoB,geoC]
scatters = [[geoA,geoB],[geoB,geoC],[geoC,geoA]]
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_alpha/'
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import A01Unzip as A01
import A02CreateCsv as A02
import A03CreateHtml as A03

dfs = []
for name in names:

    if unzip:
        print(log.getTime(), 'Unzip', name)
        A01.unzipProteome(name, dir)

    if csv:
        print(log.getTime(), 'Create CSV', name)
        A02.createCsvCorrelations(name, dir, geos, "Contacts", chunk, overwrite)

    print(log.getTime(),'load csv', "Csv/" + name + "_Contacts.csv")
    df_contacts = pd.read_csv("Csv/" + name + "_Contacts.csv")
    dfs.append([name, df_contacts])

    if html:
        print(log.getTime(), 'Create HTML', name)
        A03.createHtmlFromGeos(name, 'Contacts',df_contacts,geos,scatters)





