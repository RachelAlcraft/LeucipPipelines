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
#runs.append(['High',1000])
#runs.append(['Redo',1000])
runs.append(['Redo_GLY_IDEAL',1000])
runs.append(['High_GLY_IDEAL',1000])
#runs.append(['Redo_GLY',1000])
#runs.append(['High_GLY',1000])

#runs.append(['SYN_GLY',1000])

for tag,iters in runs:
    if run_out_proc:
        print('Starting subprocess',tag,iters)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A11_ProteinTop_OutProc.py'
        commands2 = str(tag)
        commands2 += ' ' + str(iters)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,tag,str(iters)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',tag,iters)
        import A11_ProteinTop_OutProc as inproca
        inproca.proteinTop(tag,str(iters))
