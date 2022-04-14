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
    csv_correlations = "Csv/10_DivCorr_" + tag + ".csv"
    print(csv_correlations)
    html_filename = 'Html/B01_ProteinTop_' + tag + str_iters + stat + '.html'
    log_file = "Log/B01_ProteinTop_" + tag + str_iters + stat+ ".log"
    title = 'Alcraft-Williams Associations: Most and Least Associated ' + tag
    openLog(log_file, str_iters)
    ###############################################################################################
    csvA = csv_final
    csvB = ''
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

    iters = int(str_iters)
    dataA = pd.read_csv(csvA)
    dataB = None
    if csvB != '':
        dataB = pd.read_csv(csvB)
    geos = glob.getGeos(False,True)

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        dataA = glob.trimDihs(dataA, 15)
        if csvB != '':
            dataB = glob.trimDihs(dataB, 15)

    log(log_file,'Create Williams Coefficient Maker')
    aw = None
    if csvB != '':
        aw = awa.AlcraftWilliamsAssociation(dataA,dataB, piters=iters,method=method)
    else:
        aw = awa.AlcraftWilliamsAssociation(dataA,piters=iters,method=method)
    complete = pd.read_csv(csv_correlations)
    complete = complete.sort_values(by=stat, ascending=False)
    recreate_html = True
    dont_use = ['C-1:CB-1','CA-1:C-1','N-1:C-1','N-1:CA-1:C-1','CB:CA:C','N-1:CA-1','N-1:CB-1','CA-1:CB-1']

    if recreate_html:
        used_geos = []
        log(log_file,'### Creating html reports')
        rep_mak = re.ReportExport(title, html_filename, cols=6)
        rep_mak.addLineComment('Most correlated geos')
        i = -1
        while len(used_geos) < num_top:
            i += 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos and geoA not in dont_use and geoB not in dont_use:
                log(log_file,str(len(used_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB)
                used_geos.append(geoA + '_' + geoB)
                div = aw.addAssociation([geoA,geoB])
                if csvB == '':
                    dataB = aw.getShuffledData(dataA,[geoA,geoB])

                stat, pvalue, A, D, B = div.metric, div.pvalue, div.matA, div.matDiff, div.matB
                histA,histB = div.phistA,div.phistB
                maxV = max(np.max(D), -1*np.min(B))
                rep_mak.addPlot2d(dataA, 'seaborn', title='Observed stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), geo_x=geoA, geo_y=geoB, hue='aa',palette='tab20')
                rep_mak.addPlot2d(dataB, 'scatter', title='Compare', geo_x=geoA, geo_y=geoB, hue=geoA)
                if len(histA) > 0:
                    rep_mak.addPlot1d(histA, 'histogram', title='', overlay=True, alpha=0.5,palette='steelblue')
                    rep_mak.addPlot1d(histB, 'histogram', title='', alpha=0.5,palette='Firebrick')
                else:
                    rep_mak.addBoxComment(('No histogram calculated'))
                rep_mak.addSurface(A, 'Original Data', palette='Blues', colourbar=False)
                rep_mak.addSurface(D, 'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)),  cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                rep_mak.addSurface(B, 'Compare Data', palette='Reds', colourbar=False)

        log(log_file,'############# Least correlated ###########################')
        rep_mak.addLineComment('Least correlated geos')
        i = len(complete.index)
        last_geos = []
        while len(last_geos) < num_bottom:
            i -= 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos and geoA not in dont_use and geoB not in dont_use:
                log(log_file, str(len(used_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB)
                last_geos.append(geoA + '_' + geoB)
                div = aw.addAssociation([geoA, geoB])
                if csvB == '':
                    dataB = aw.getShuffledData(dataA, [geoA, geoB])

                stat, pvalue, A, D, B = div.metric, div.pvalue, div.matA, div.matDiff, div.matB
                histA, histB = div.phistA, div.phistB
                maxV = max(np.max(D), -1*np.min(B))
                rep_mak.addPlot2d(dataA, 'seaborn',title='Observed stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), geo_x=geoA, geo_y=geoB, hue='aa', palette='tab20')
                rep_mak.addPlot2d(dataB, 'scatter', title='Compare', geo_x=geoA, geo_y=geoB, hue=geoA)
                if len(histA) > 0:
                    rep_mak.addPlot1d(histA, 'histogram', title='', overlay=True, alpha=0.5, palette='steelblue')
                    rep_mak.addPlot1d(histB, 'histogram', title='', alpha=0.5, palette='Firebrick')
                else:
                    rep_mak.addBoxComment(('No histogram calculated'))
                rep_mak.addSurface(A, 'Original Data', palette='Blues', colourbar=False)
                rep_mak.addSurface(D,'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                rep_mak.addSurface(B, 'Compare Data', palette='Reds', colourbar=False)

        log(log_file, 'Finally print out to ' + html_filename)
        rep_mak.printReport()

###################################################################################
if __name__ == '__main__':
    globals()['proteinTop'](sys.argv[1],sys.argv[2],sys.argv[3])




