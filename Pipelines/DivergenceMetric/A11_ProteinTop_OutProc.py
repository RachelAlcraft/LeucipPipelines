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
from LeucipPy import WilliamsDivergenceMaker as wcm

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
def proteinTop(tag,str_iters):
    csv_final = "Csv/PW_" + tag + "_01_Geometry.csv"
    tag_corr = tag
    if 'IDEAL' in tag:
        tag_corr = tag[:len(tag)-6]
    csv_correlations = "Csv/10_DivCorr_" + tag_corr + ".csv"
    print(csv_correlations)
    html_filename = 'Html/11_ProteinTop_' + tag + str_iters + '.html'
    log_file = "Log/11_ProteinTop_" + tag + str_iters + ".log"
    title = 'Williams Divergence from Trivial: Most and Least Correlated ' + tag
    openLog(log_file, str_iters)
    ###############################################################################################
    iters = int(str_iters)
    data = pd.read_csv(csv_final)
    geos = []
    for col in data.columns:
        #if ':' in col and 'info' not in col and '(' not in col and '{' not in col:
        #if ':' in col and 'info' not in col and '{' not in col:
        if ':' in col and 'info' not in col:
            geos.append(col)

    #geos = ['N:O','N:CA','CA:C','N:CA:C:N+1']

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        data = glob.trimDihs(data, 15)

    log(log_file,'Create Williams Coefficient Maker')
    wcc = wcm.WilliamsDivergenceMaker(data,geos,density=density,log=1,norm=True,pval_iters=iters,delay_load=True)

    complete = pd.read_csv(csv_correlations)
    complete = complete.sort_values(by='stat', ascending=False)
    recreate_html = True
    if recreate_html:
        used_geos = []
        log(log_file,'### Creating html reports')
        rep_mak = hrm.HtmlReportMaker(title, html_filename, cols=6)
        rep_mak.addLineComment('Most correlated geos')
        i = -1
        while len(used_geos) < num_top:
            i += 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos:
                log(log_file,str(len(used_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB)
                used_geos.append(geoA + '_' + geoB)
                div = wcc.getCorrelation([geoA,geoB])
                cm_data = wcc.data
                cm_data = cm_data.sort_values(by='aa+1')

                data_ab = glob.trimGeos(cm_data, 15, geoA, geoB)
                df_rand = wcc.randomiseData(cm_data[[geoA, geoB]])
                df_rand = glob.trimGeos(df_rand, 15, geoA, geoB)

                stat, pvalue, A, D, B = div.stat, div.p_value, div.histAB, div.diffAB, div.convAB
                hist = div.p_hist
                maxV = max(np.max(A), np.max(B))
                rep_mak.addPlot2d(data_ab, 'seaborn', title='Observed Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), geo_x=geoA, geo_y=geoB, hue='aa',palette='tab20')
                rep_mak.addPlot2d(df_rand, 'scatter', title='rand', geo_x=geoA, geo_y=geoB, hue=geoA)
                if len(hist['divergence_shuffled']) > 0:
                    rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence_shuffled', title='', overlay=True, alpha=0.5,palette='steelblue')
                    rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence_resampled', title='', alpha=0.5,palette='Firebrick')
                else:
                    rep_mak.addBoxComment(('No histogram calculated'))
                rep_mak.addSurface(A, 'Original Data', cmin=0, cmax=maxV, palette='Blues', colourbar=False)
                rep_mak.addSurface(D, 'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)),  cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                rep_mak.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds', colourbar=False)

        log(log_file,'############# Least correlated ###########################')
        rep_mak.addLineComment('Least correlated geos')
        i = len(complete.index)
        last_geos = []
        while len(last_geos) < num_bottom:
            i -= 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]

            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos:
                log(log_file,str(len(last_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB )
                last_geos.append(geoA + '_' + geoB)
                div = wcc.getCorrelation([geoA, geoB])
                cm_data = wcc.data
                cm_data = cm_data.sort_values(by='aa+1')

                data_ab = glob.trimGeos(cm_data, 15, geoA, geoB)
                df_rand = wcc.randomiseData(cm_data[[geoA, geoB]])
                df_rand = glob.trimGeos(df_rand, 15, geoA, geoB)

                stat, pvalue, A, D, B = div.stat, div.p_value, div.histAB, div.diffAB, div.convAB
                hist = div.p_hist
                maxV = max(np.max(A), np.max(B))
                rep_mak.addPlot2d(data_ab, 'seaborn', title='Observed Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), geo_x=geoA, geo_y=geoB, hue='aa',palette='tab20_r')
                rep_mak.addPlot2d(df_rand, 'scatter', title='rand', geo_x=geoA, geo_y=geoB, hue=geoA)
                if len(hist['divergence_shuffled']) > 0:
                    rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence_shuffled', title='', overlay=True, alpha=0.5,palette='steelblue')
                    rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence_resampled', title='', alpha=0.5,palette='Firebrick')
                else:
                    rep_mak.addBoxComment(('No histogram calculated'))
                rep_mak.addSurface(A, 'Original Data', cmin=0, cmax=maxV, palette='Blues', colourbar=False)
                rep_mak.addSurface(D, 'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                rep_mak.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds', colourbar=False)

        log(log_file, 'Finally print out to ' + html_filename)
        rep_mak.printReport()

###################################################################################
if __name__ == '__main__':
    globals()['proteinTop'](sys.argv[1],sys.argv[2])




