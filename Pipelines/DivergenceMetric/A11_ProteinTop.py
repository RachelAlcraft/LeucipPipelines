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
#runs.append(['Redo_GLY',0,True,True,False])
#runs.append(['Redo_GLY',1000,True,False,True])
#runs.append(['High_GLY',0,True,True,False])
#runs.append(['High_GLY',1000,True,False,True])
#runs.append(['Redo',0,True,True,False])
#runs.append(['Redo',1000,True,False,True])
runs.append(['High',0,True,True,False])
runs.append(['High',1000,True,False,True])

for tag,iters,modify_csv,recreate_divergence,recreate_html in runs:
    if run_out_proc:
        print('Starting subprocess',tag,iters,modify_csv,recreate_divergence,recreate_html)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A11_ProteinTop_OutProc.py'
        commands2 = str(tag)
        commands2 += ' ' + str(iters)
        commands2 += ' ' + str(modify_csv)
        commands2 += ' ' + str(recreate_divergence)
        commands2 += ' ' + str(recreate_html)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,tag,str(iters),str(modify_csv),str(recreate_divergence),str(recreate_html)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',tag,iters,modify_csv,recreate_divergence,recreate_html)
        import A11_ProteinTop_OutProc as inproca
        inproca.proteinTop(tag,str(iters),str(modify_csv),str(recreate_divergence),str(recreate_html))
