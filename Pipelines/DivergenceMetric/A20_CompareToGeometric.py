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
#runs.append(['High_GLY','High_GLY_IDEAL','',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL','',False,True])
#runs.append(['Redo_GLY','Redo_GLY_IDEAL','',False,False])
#runs.append(['Redo_GLY','Redo_GLY_IDEAL','',False,True])
#runs.append(['High_GLY','High_GLY','High_GLY_IDEAL_PAIRED','',False,0])
#runs.append(['High_GLY','High_GLY','High_GLY_IDEAL_UNPAIRED','',False,0])
#runs.append(['High_GLY','High_GLY_IDEAL_PAIRED','High_GLY_IDEAL_UNPAIRED','',False,3])
#runs.append(['High_GLY','High_GLY','High_GLY_IDEAL_PAIRED','',False,'stat'])
#runs.append(['High','High','High_IDEAL_PAIRED','',False,'stat'])

runs.append(['High_GLY','High_GLY','High_GLY_IDEAL_PAIRED','',False,'stat'])
runs.append(['High','High','High_IDEAL_PAIRED','',False,'stat'])
runs.append(['High_GLY','High_GLY','High_GLY_IDEAL_PAIRED','N:CA:C',False,'stat'])
runs.append(['High','High','High_IDEAL_PAIRED','N:CA:C',False,'stat'])
runs.append(['High_GLY','High_GLY','High_GLY_IDEAL_PAIRED','C:O',False,'stat'])
runs.append(['High','High','High_IDEAL_PAIRED','C:O',False,'stat'])



#runs.append(['High_GLY','High_GLY_IDEAL_PAIRED','N:CA:C',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_UNPAIRED','N:CA:C',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_SAMPLED','N:CA:C',False,False])

#runs.append(['High_GLY','High_GLY_IDEAL_SAMPLED','C:O',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_PAIRED','C:O',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_PAIRED','CA:C:N+1',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_UNPAIRED','C:O',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_UNPAIRED','CA:C:N+1',False,False])
#runs.append(['High_GLY','High_GLY_IDEAL_SAMPLED','CA:C:N+1',False,False])



for ID,tag1,tag2, geoInc,dssp,sort_geo in runs:
    if run_out_proc:
        print('Starting subprocess',ID,tag1,tag2,dssp,sort_geo)
        exe = sys.executable
        command = 'C:/Dev/Github/LeucipPipelines/Pipelines/DivergenceMetric/A20_CompareToGeometric_OutProc.py'
        commands2 = str(ID)
        commands2 += str(tag1)
        commands2 += ' ' + str(tag2)
        commands2 += ' ' + str(geoInc)
        commands2 += ' ' + str(dssp)
        commands2 += ' ' + str(sort_geo)
        print('"' + exe + '"' ,command,commands2)
        pigP = sub.Popen([exe,command,ID,tag1,tag2,geoInc,str(dssp),str(sort_geo)], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
    else:
        print('Running in process',ID,tag1,tag2,geoInc,dssp,sort_geo)
        import A20_CompareToGeometric_OutProc as inproca
        inproca.proteinTopCompare(ID,tag1,tag2,geoInc,str(dssp),str(sort_geo))
