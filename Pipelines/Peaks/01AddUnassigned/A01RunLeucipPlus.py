
import os
import subprocess as sub
import shutil

def runCreateElectronDensity(pdbcodes,ExePath,LogDirectory):
    print('LeucipPipeline:01Images - 2) Creating electron density -----------------------------------------------')
    count = 0
    for pdbcode in pdbcodes:
        count += 1
        print("LeuciPipelines:Embellish: running LeucipPlus.exe", count, '/', len(pdbcodes),pdbcode)
        cmd = "EMBELLISHFILE|" + pdbcode + "|5|1|0|GRIDSIZE|CENTRAL|LINEAR|PLANAR|" + LogDirectory
        print(cmd)
        #movepath = 'Output/' + set + '_' + id + '_' + pdb + '_SLICESFILE.csv'
        #if not os.path.exists(movepath):
        pigP = sub.Popen([ExePath, cmd], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
        #outputpath = LogDirectory + pdb + '_SLICESFILE.csv'
        #shutil.move(outputpath, movepath)