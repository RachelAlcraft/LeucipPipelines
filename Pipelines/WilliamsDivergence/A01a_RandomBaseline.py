import sys
import subprocess as sub

run_out_proc = True

runs2 = []
#runs2.append([500,0.5])
runs2.append([0,1,1])
runs2.append([0,1,10])
runs2.append([0,1,20])
runs2.append([0,1,100])
#runs2.append([500,5])
#runs2.append([500,10])
#runs2.append([500,20])

for sample_size,density,vari in runs2:
    if run_out_proc:
        print('Starting subprocess',sample_size,density)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/WilliamsDivergence/A01a_RandomBaseline_OutProc.py'
        commands2 = str(sample_size)
        commands2 += ' ' + str(density)
        commands2 += ' ' + str(vari)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,str(sample_size),str(density),str(vari)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process', sample_size,density)
        import A01a_RandomBaseline_OutProc as inproca
        inproca.randomBaselineOne(str(sample_size),str(density),str(vari))
