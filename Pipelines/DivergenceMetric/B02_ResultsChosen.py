import sys
import subprocess as sub

run_out_proc = True

# INPUTS #
#control which steps you want to run
'''

'''
# stat_trivial
# stat_obs_paired
# stat_obs_unpaired
# stat_paired_unpaired
# kl_trivial
# kl_obs_paired
# kl_obs_unpaired
# kl_paired_unpaired

runs = []
#runs.append(['High',100,'diff'])
#runs.append(['High',10,'stat_obs_unpaired'])
runs.append(['High',10,'stat_obs_paired'])
#runs.append(['High',10,'stat_paired_unpaired'])
#runs.append(['High',10,'kl_trivial'])
#runs.append(['High',10,'kl_obs_unpaired'])
#runs.append(['High',10,'kl_paired_unpaired'])

for tag,iters,stat in runs:
    if run_out_proc:
        print('Starting subprocess',tag,iters)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/B02_ResultsChosen_OutProc.py'
        commands2 = str(tag)
        commands2 += ' ' + str(iters)
        commands2 += ' ' + str(stat)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,tag,str(iters),stat], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',tag,iters,stat)
        import B02_ResultsChosen_OutProc as inproca
        inproca.proteinTop(tag,str(iters),stat)
