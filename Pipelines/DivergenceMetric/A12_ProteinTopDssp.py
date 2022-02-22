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
#runs.append('Redo_GLY')
#runs.append('Redo')
runs.append('High_GLY')
#runs.append('High')
#runs.append('SYN_GLY')

for tag in runs:
    if run_out_proc:
        print('Starting subprocess',tag)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A12_ProteinTopDssp_OutProc.py'
        commands2 = str(tag)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,tag], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',tag)
        import A12_ProteinTopDssp_OutProc as inproca
        inproca.proteinTopDssp(tag)
