

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
def runCreateIAM(pdbName,ExePath,LogDirectory):
    print('LeucipPipeline: -Creating IAM density -----------------------------------------------')
    cmd = 'DEFORMATIONFILE|' + pdbName + '|5|2|-1|x|x|' + LogDirectory + '|' + LogDirectory + '|' +LogDirectory
    print('command=',cmd)
    pigP = sub.Popen([ExePath, cmd], stdout=sub.PIPE)
    resultP = pigP.communicate(input=b"This is sample text.\n")
    exe_resultP = str(resultP[0], 'utf-8')
    pigP.kill()





#***************************************************************************************************************
