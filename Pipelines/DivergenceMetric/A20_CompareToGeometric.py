import sys
import subprocess as sub

run_out_proc = True

# INPUTS #
#control which steps you want to run
'''
modify_csv,recreate_divergence,recreate_html = True,False,True
tag = 'Redo' #Redo or High or Redo_GLY or High_GlY
iters = 1000
'''
runs = []
runs.append(['High_GLY','SYN_GLY','N:CA:C'])
runs.append(['High_GLY','SYN_GLY',''])
runs.append(['Redo_GLY','SYN_GLY','N:CA:C'])
runs.append(['Redo_GLY','SYN_GLY',''])
runs.append(['Redo_GLY','High_GLY',''])


for tag1,tag2, geoInc in runs:
    if run_out_proc:
        print('Starting subprocess',tag1,tag2)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A20_CompareToGeometric_OutProc.py'
        commands2 = str(tag1)
        commands2 += ' ' + str(tag2)
        commands2 += ' ' + str(geoInc)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,tag1,tag2,geoInc], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',tag1,tag2,geoInc)
        import A20_CompareToGeometric_OutProc as inproca
        inproca.proteinTopCompare(tag1,tag2,geoInc)
