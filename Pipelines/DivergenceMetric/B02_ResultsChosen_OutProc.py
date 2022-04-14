'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
density = 5
num_top = 40
num_bottom = 10
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help
import A0_Globals as glob

from LeucipPy import HtmlReportMaker as hrm
from nDimAssociations import AlcraftWilliamsAssociation as awa
from nDimAssociations import ReportExport as re

import pandas as pd
import numpy as np
from datetime import datetime as dt

########################################################################################
def getTime():
    now = dt.now()
    tm =  now.strftime("%a %d%h%y-%H:%M:%S")
    return tm
def openLog(logfile, msg):
    f = open(logfile, "w")
    f.write(getTime() + " :\tStarting log file for LeucipPy\n" + msg + '\n')
    f.close()
    print(msg)
def log(logfile, msg):
    f = open(logfile, "a")
    f.write(getTime() + " :\t" + msg + '\n')
    f.close()
    print(msg)
########################################################################################
def proteinTop(tag,str_iters,stat):
    csv_final = "Csv/PW_" + tag + "_01_Geometry.csv"
    csv_ideal_paired = "Csv/PW_" + tag + "_IDEAL_PAIRED_01_Geometry.csv"
    csv_ideal_unpaired = "Csv/PW_" + tag + "_IDEAL_UNPAIRED_01_Geometry.csv"
    html_filename = 'Html/B02_ProteinChosen_' + tag + str_iters + stat + '.html'
    log_file = "Log/B02_Chosen_" + tag + str_iters + stat+ ".log"
    title = 'Alcraft-Williams Associations: Chosen Associations ' + tag
    openLog(log_file, str_iters)
    ###############################################################################################
    csvA = csv_final
    csvP = csv_ideal_paired
    csvU = csv_ideal_unpaired
    method = stat
    '''
    method = 'diff'
    if stat == 'stat_obs_paired':
        csvB = csv_ideal_paired
    elif stat == 'stat_obs_unpaired':
        csvB = csv_ideal_unpaired
    elif stat == 'stat_paired_unpaired':
        csvA = csv_ideal_paired
        csvB = csv_ideal_unpaired
    elif stat == 'kl_trivial':
        method = 'k-l'
    elif stat == 'kl_obs_paired':
        csvB = csv_ideal_paired
        method = 'k-l'
    elif stat == 'kl_obs_unpaired':
        csvB = csv_ideal_unpaired
        method = 'k-l'
    elif stat == 'kl_paired_unpaired':
        csvA = csv_ideal_paired
        csvB = csv_ideal_unpaired
        method = 'k-l'
    '''

    iters = int(str_iters)
    dataA = pd.read_csv(csvA)
    dataP = pd.read_csv(csvP)
    dataU = pd.read_csv(csvU)

    geos = glob.getGeos(False,True)

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        dataA = glob.trimDihs(dataA, 15)
        dataP = glob.trimDihs(dataP, 15)
        dataU = glob.trimDihs(dataU, 15)

    log(log_file,'Create Williams Coefficient Maker')
    awA = awa.AlcraftWilliamsAssociation(dataA,piters=iters,method=method)
    awP = awa.AlcraftWilliamsAssociation(dataA,dataP,piters=iters,method=method)
    awU = awa.AlcraftWilliamsAssociation(dataA,dataU, piters=iters, method=method)
    awPU = awa.AlcraftWilliamsAssociation(dataP,dataU, piters=iters, method=method)
    awPr = awa.AlcraftWilliamsAssociation(dataP, piters=iters, method=method)
    awUr = awa.AlcraftWilliamsAssociation(dataU, piters=iters, method=method)


    recreate_html = True
    chosen = [  ['N:O','O:CB'],
                ['N-1:C-1', 'C-1:CB-1'],
                ['N+1:O+1','O+1:CB+1'],
                ['N-1:O-1', 'O-1:CB-1'],
                ['C-1:N:CA:C','N:CA:C:N+1'],
                ['CA:C:O','O:C:N+1'],
                ['CA:C:O','CA:C:N+1'],
                ['O:C:N+1', 'CA:C:N+1'],
                ['CA-1:N+1', 'CA-1:CA+1', 'CA-1:CA:CA+1']
                ]

    if recreate_html:
        log(log_file,'### Creating html reports')
        rep_mak = re.ReportExport(title, html_filename, cols=6)
        for geos in chosen:
            if True:
                log(log_file,str(geos) + '.........')
                rep_mak.addLineComment(str(geos))
                divA = awA.addAssociation(geos)
                divP = awP.addAssociation(geos)
                divU = awU.addAssociation(geos)
                divPr = awPr.addAssociation(geos)
                divUr = awUr.addAssociation(geos)
                divPU = awPU.addAssociation(geos)

                dataAs = awA.getShuffledData(dataA, geos)
                dataPs = awA.getShuffledData(dataP, geos)
                dataUs = awA.getShuffledData(dataU, geos)

                statA, pvalueA, Aa, Da, Ba = divA.metric, divA.pvalue, divA.matA, divA.matDiff, divA.matB
                statP, pvalueP, Ap, Dp, Bp = divP.metric, divP.pvalue, divP.matA, divP.matDiff, divP.matB
                statU, pvalueU, Au, Du, Bu = divU.metric, divU.pvalue, divU.matA, divU.matDiff, divU.matB
                statPr, pvaluePr, Apr, Dpr, Bpr = divPr.metric, divPr.pvalue, divPr.matA, divPr.matDiff, divPr.matB
                statUr, pvalueUr, Aur, Dur, Bur = divUr.metric, divUr.pvalue, divUr.matA, divUr.matDiff, divUr.matB
                statPU, pvaluePU, Apu, Dpu, Bpu = divPU.metric, divPU.pvalue, divPU.matA, divPU.matDiff, divPU.matB
                histA,histB = divPU.phistA,divPU.phistB
                maxV = max(np.max(Dpu), -1*np.min(Dpu))
                geoA = geos[0]
                geoB = geos[1]
                hue = geoA
                if len(geos)>2:
                    hue = geos[2]
                rep_mak.addLineComment('Un-Paired  stat=' + str(round(statPU, 3)) + ' pvalue=' + str(round(pvaluePU, 3)))
                rep_mak.addPlot2d(dataA, 'scatter', title='Observed stat=' + str(round(statA, 3)) + ' pvalue=' + str(round(pvalueA, 3)), geo_x=geoA, geo_y=geoB, hue=hue,palette='Spectral')
                rep_mak.addPlot2d(dataAs, 'scatter', title='Randomised', geo_x=geoA, geo_y=geoB, hue=geoA)
                rep_mak.addPlot2d(dataP, 'scatter', title='Obs-Paired stat=' + str(round(statP, 3)) + ' pvalue=' + str(round(pvalueP, 3)), geo_x=geoA, geo_y=geoB, hue=hue, palette='Spectral')
                rep_mak.addPlot2d(dataPs, 'scatter', title='Paired stat=' + str(round(statPr, 3)) + ' pvalue=' + str(round(pvaluePr, 3)), geo_x=geoA, geo_y=geoB, hue=geoA)
                rep_mak.addPlot2d(dataU, 'scatter', title='Obs-Unpaired stat=' + str(round(statU, 3)) + ' pvalue=' + str(round(pvalueU, 3)), geo_x=geoA, geo_y=geoB, hue=hue, palette='Spectral')
                rep_mak.addPlot2d(dataUs, 'scatter', title='Unpaired stat=' + str(round(statUr, 3)) + ' pvalue=' + str(round(pvalueUr, 3)), geo_x=geoA, geo_y=geoB, hue=geoA)
                #if len(histA) > 0:
                #    rep_mak.addPlot1d(histA, 'histogram', title='', overlay=True, alpha=0.5,palette='steelblue')
                #    rep_mak.addPlot1d(histB, 'histogram', title='Un-Paired stat=' + str(round(statPU, 3)) + ' pvalue=' + str(round(pvaluePU, 3)), alpha=0.5,palette='Firebrick')
                #else:
                #    rep_mak.addBoxComment(('No histogram calculated'))
                #rep_mak.addSurface(A, 'Original Data', palette='Blues', colourbar=False)
                #rep_mak.addSurface(D, 'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)),  cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                #rep_mak.addSurface(B, 'Convolved Data', palette='Reds', colourbar=False)

        log(log_file,'############# Least correlated ###########################')
        rep_mak.addLineComment('Least correlated geos')

        log(log_file, 'Finally print out to ' + html_filename)
        rep_mak.printReport()

###################################################################################
if __name__ == '__main__':
    globals()['proteinTop'](sys.argv[1],sys.argv[2],sys.argv[3])




