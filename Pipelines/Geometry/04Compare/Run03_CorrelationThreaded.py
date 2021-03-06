#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import subprocess as sub

####  User inputs  #########################################################
order_geo_list,make_html = True,True
runs = [] #ID,csv,geoA,aa_inc,aa_exc,hue,tag,over_geos,chunk,outlir_cut

#runs.append(['Correlation_OCAm1_','PW_High_GLY_02_Geometry.csv','O:CA-1',['GLY'],'','N:CA:C','',[],45,25])
#runs.append(['Correlation_OCAm1_','PW_High_02_Geometry.csv','O:CA-1',[''],'PRO','N:CA:C','',[],45,25])

#runs.append(['Correlation_NO_','PW_High_02_Geometry.csv','N:O',[''],'PRO','bfactor','',[],45,25])
#runs.append(['Correlation_NO_','PW_High_02_Geometry.csv','N:O',['PRO','ASP','LYS','TYR'],'','bfactor','',[],45,25])
#runs.append(['Correlation_NO_','PW_High_GLY_02_Geometry.csv','N:O',['GLY'],'','bfactor','',[],45,25])

#runs.append(['Correlation_PHI_','PW_High_02_Geometry.csv','C-1:N:CA:C',[''],'PRO','bfactor','',[],45,25])
#runs.append(['Correlation_PHI_','PW_High_02_Geometry.csv','C-1:N:CA:C',['PRO','ASP','LYS','TYR'],'','bfactor','',[],45,25])
#runs.append(['Correlation_PHI_','PW_High_GLY_02_Geometry.csv','C-1:N:CA:C',['GLY'],'','O:CA-1','',[],45,25])

runs.append(['Correlation_PSI_','PW_High_02_Geometry.csv','N:CA:C:N+1',[''],'PRO','bfactor','',[],45,25])
runs.append(['Correlation_PSI_','PW_High_02_Geometry.csv','N:CA:C:N+1',['PRO','ASP','LYS','TYR'],'','bfactor','',[],45,25])
runs.append(['Correlation_PSI_','PW_High_GLY_02_Geometry.csv','N:CA:C:N+1',['GLY'],'','bfactor','',[],45,25])

#runs.append(['Correlation_TAU_','PW_High_02_Geometry.csv','N:CA:C',[''],'PRO','bfactor','',[],45,25])
#runs.append(['Correlation_TAU_','PW_High_02_Geometry.csv','N:CA:C',['PRO','ASP','LYS','TYR'],'','bfactor','',[],45,25])
#runs.append(['Correlation_TAU_','PW_High_GLY_02_Geometry.csv','N:CA:C',['GLY'],'','bfactor','',[],45,25])

#runs.append(['Correlation_TAUm1_','PW_High_02_Geometry.csv','C-1:N:CA',[''],'PRO','bfactor','',[],45,25])
#runs.append(['Correlation_TAUm1_','PW_High_02_Geometry.csv','C-1:N:CA',['PRO','ASP','LYS','TYR'],'','bfactor','',[],45,25])
#runs.append(['Correlation_TAUm1_','PW_High_GLY_02_Geometry.csv','C-1:N:CA',['GLY'],'','bfactor','',[],45,25])

#runs.append(['Correlation_TAUp1_','PW_High_02_Geometry.csv','CA:C:N+1',[''],'PRO','bfactor','',[],45,25])
#runs.append(['Correlation_TAUp1_','PW_High_02_Geometry.csv','CA:C:N+1',['PRO','ASP','LYS','TYR'],'','bfactor','',[],45,25])
#runs.append(['Correlation_TAUp1_','PW_High_GLY_02_Geometry.csv','CA:C:N+1',['GLY'],'','bfactor','',[],45,25])

geos_to_abs = ['CA:C:O:N+1','CA-1:C-1:N:CA','CA:C:N+1:CA+1']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Some filtering of the data
aa_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
aa_NO_gly = ['ALA','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###   THREAD WORK    #########################################
### https://www.tutorialspoint.com/python/python_multithreading.htm
### https://www.toptal.com/python/beginners-guide-to-concurrency-and-parallelism-in-python


##  MAIN   ###################################################

import pandas as pd
import sys
from datetime import datetime
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Helpers as hlp
import A01MakeCsv as A01
import A02CreateHtml as A02
import A03FilterAndRandomise as A03
import A04DifferenceImage as A04
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
startA = datetime.now()
startB = datetime.now()

run_count = 0
for ID,csv,geoA,aa_inc,aa_exc,hue,tag,over_geos,chunk,outlier_cut in runs:
    run_count += 1
    csv_filename = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/" + csv
    print(run_count,'/',len(runs),ID,geoA,aa_inc,aa_exc,csv)
    df_geometryAllAA = pd.read_csv(csv_filename)
    for aa_inc in aa_inc:
        if aa_inc != '':
            df_geometryAll = df_geometryAllAA.query("aa == '" + aa_inc + "'")
        if aa_exc != '':
            df_geometryAll = df_geometryAllAA.query("aa != '" + aa_exc + "'")

        for gabs in geos_to_abs:
            df_geometryAll[gabs] =abs(df_geometryAll[gabs])

        #cut off the extreme outliers
        df_geometryAll = df_geometryAll.sort_values(by=geoA, ascending=True)
        df_geometryAll = df_geometryAll.iloc[outlier_cut:,:]
        df_geometryAll = df_geometryAll.sort_values(by=geoA, ascending=False)
        df_geometryAll = df_geometryAll.iloc[outlier_cut:,:]

        cols = []
        for col in df_geometryAll.columns:
            if ':' in col and 'info' not in col:
                cols.append(col)

        df_geometryA = df_geometryAll[cols]
        cols.append(hue)
        df_geometryB = df_geometryAll[cols]

        if order_geo_list:
            geo_list = df_geometryA.columns
            ordered_list_by_stats = []
            count = 0
            for geoB in geo_list:
                count += 1
                if geoA != geoB:
                    # First make the appropriate distributions
                    print('Make distribution difference for', geoB, count, '/', len(df_geometryA.columns))
                    df_col = df_geometryA.sort_values(by=geoB, ascending=True)
                    df_col = df_col.iloc[outlier_cut:, :]
                    df_col = df_col.sort_values(by=geoB, ascending=False)
                    df_col = df_col.iloc[outlier_cut:, :]
                    arA, arB, arDiff, minv, maxv, stat,arDiffSq = A04.createIdealisedDifferencePlot(geoA, geoB, df_col)
                    ordered_list_by_stats.append([geoB, stat])

            print('sort list to most dependent at the top')
            ordered_list_by_stats = sorted(ordered_list_by_stats, key=lambda tup: tup[1], reverse=True)

            ordered_dic = {}
            ordered_dic['geo'] = []
            ordered_dic['stat'] = []
            for geoB, stat in ordered_list_by_stats:
                ordered_dic['geo'].append(geoB)
                ordered_dic['stat'].append(stat)

            df_ordered = pd.DataFrame.from_dict(ordered_dic)
            df_ordered.to_csv("C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/OrderedGeos_" + ID + aa_inc + ".csv", index=False)

        if make_html:
            ordered_geos_csv = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/OrderedGeos_" + ID + aa_inc + ".csv"

            df_ordered = pd.read_csv(ordered_geos_csv)
            allgeos = df_ordered['geo'].values

            start_count, end_count = 0, 0 + chunk
            while start_count < len(allgeos):
                if end_count >  len(allgeos):
                    end_count = len(allgeos)
                html_filename = 'C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/ThreadHtml/' + ID + aa_inc + tag + '_Dependency_' + str(start_count) + '_' + str(end_count) + '.html'

                    # https://stackoverflow.com/questions/2046603/is-it-possible-to-run-function-in-a-subprocess-without-threading-or-writing-a-se
                print('Starting subprocess', geoA, start_count, end_count, len(allgeos))
                #exe = "C:/Program Files (x86)/Microsoft Visual Studio/Shared/Python37_64/python.exe"
                exe = sys.executable
                command = 'C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/ExtHtmlReport.py'
                #commands2 = ' runOneHtml'
                commands2 = ordered_geos_csv
                commands2 += ' ' + csv_filename
                commands2 += ' ' + html_filename
                commands2 += ' ' + geoA
                commands2 += ' ' + str(start_count)
                commands2 += ' ' + str(end_count)
                commands2 += ' ' + aa_inc
                commands2 += ' ' + aa_exc
                commands2 += ' ' + hue
                print('"' + exe + '"' ,command,commands2)
                pigP = sub.Popen([exe,command,ordered_geos_csv,csv_filename,html_filename,geoA,str(start_count),str(end_count),aa_inc,aa_exc,hue], stdout=sub.PIPE)
                resultP = pigP.communicate(input=b"This is sample text.\n")
                exe_resultP = str(resultP[0], 'utf-8')
                pigP.kill()

                start_count += chunk
                end_count += chunk
                end = datetime.now()
                hlp.printTime(startB, end)
                startB = datetime.now()


print('Total time was...')
end = datetime.now()
hlp.printTime(startA,end)


