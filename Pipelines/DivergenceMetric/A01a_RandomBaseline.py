import sys
import subprocess as sub

run_out_proc = False

runs2 = []
runs2.append([0,0.5,1])
#runs2.append([0,1,1])
#runs2.append([0,5,1])
#runs2.append([0,10,1])


for iters,density,vari in runs2:
    if run_out_proc:
        print('Starting subprocess',iters,density)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A01a_RandomBaseline_OutProc.py'
        commands2 = str(iters)
        commands2 += ' ' + str(density)
        commands2 += ' ' + str(vari)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,str(iters),str(density),str(vari)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process', iters,density)
        import A01a_RandomBaseline_OutProc as inproca
        inproca.randomBaselineOne(str(iters),str(density),str(vari))
