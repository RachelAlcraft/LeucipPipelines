import sys
import subprocess as sub

run_out_proc = False

# INPUTS #
#control which steps you want to run
runs = []
runs.append([0.5,5,500])
runs.append([1,5,500])
runs.append([2,5,500])
runs.append([5,5,500])
runs.append([10,5,500])
runs.append([50,5,500])
#runs.append([0,5,500])


for density,bins,iters in runs:
    if run_out_proc:
        print('Starting subprocess',density,iters)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A01_CompareBaselinesVariance_OutProc.py'
        commands2 = str(density)
        commands2 += ' ' + str(bins)
        commands2 += ' ' + str(iters)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,str(density),str(bins),str(iters)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',density,iters)
        import A01_CompareBaselinesVariance_OutProc as inproca
        inproca.run(str(density),str(bins),str(iters))
