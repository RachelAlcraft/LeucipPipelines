import sys
import subprocess as sub

runs = []
#runs.append([True,50000])
#runs.append([True,500])
#runs.append([True,5000])
runs.append([True,200])

#runs.append([False,200])




for normed_corr,sample_size in runs:
    print('Starting subprocess',normed_corr,sample_size)
    exe = sys.executable
    command = 'C:/Dev/Github/LeucipPipelines/Pipelines/WilliamsDivergence/00_SyntheticRun_OutProc.py'
    commands2 = str(normed_corr)
    commands2 += ' ' + str(sample_size)
    print('"' + exe + '"' ,command,commands2)
    pigP = sub.Popen([exe,command,str(normed_corr),str(sample_size)], stdout=sub.PIPE)
    #pigP = sub.Popen([exe, command], stdout=sub.PIPE)
    resultP = pigP.communicate(input=b"This is sample text.\n")
    exe_resultP = str(resultP[0], 'utf-8')
    pigP.kill()
