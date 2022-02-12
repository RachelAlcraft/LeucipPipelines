import sys
import subprocess as sub

run_out_proc = False

runs1 = []
#runs1.append(['rand',500])
#runs1.append(['line',500])
runs1.append(['covar',500])
#runs1.append(['triangle',500])


for style,sample_size in runs1:
    if run_out_proc:
        print('Starting subprocess',style,sample_size)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/WilliamsDivergence/A01_RandomBaseline_OutProc.py'
        commands2 = str(style)
        commands2 += ' ' + str(sample_size)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,str(style),str(sample_size)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process', style, sample_size)
        import A01_RandomBaseline_OutProc as inproc
        inproc.randomBaseline(style,str(sample_size))

