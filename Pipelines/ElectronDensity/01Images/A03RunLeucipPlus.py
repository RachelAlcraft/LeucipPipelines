'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import os
import subprocess as sub
import shutil

#***************************************************************************************************************
def runCreateElectronDensity(df_leuci,ExePath,LogDirectory):
    print('LeucipPipeline:01Images - 2) Creating electron density -----------------------------------------------')
    idxs = df_leuci.index
    count = 0
    for idx in idxs:
        count += 1
        print("LeuciPipelines:01Images: running LeucipPlus.exe", count, '/', len(idxs))
        set = df_leuci['SET'][idx].lower()
        id = df_leuci['ID'][idx].lower()
        pdb = df_leuci['pdb_code'][idx].lower()
        cmd = df_leuci['CMD'][idx]
        movepath = 'Output/' + set + '_' + id + '_' + pdb + '_SLICESFILE.csv'
        if not os.path.exists(movepath):
            pigP = sub.Popen([ExePath, cmd], stdout=sub.PIPE)
            resultP = pigP.communicate(input=b"This is sample text.\n")
            exe_resultP = str(resultP[0], 'utf-8')
            pigP.kill()
            outputpath = LogDirectory + pdb + '_SLICESFILE.csv'
            shutil.move(outputpath, movepath)




#***************************************************************************************************************
