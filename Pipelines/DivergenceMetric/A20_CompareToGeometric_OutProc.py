'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
density = 5
num_top = 40
num_bottom = 10
iters = 10
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help

from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsDivergenceMaker as wcm

import pandas as pd
import numpy as np
import seaborn as sns
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
def proteinTopCompare(tagA,tagB,geoInc,str_dssp,str_sort_geo):
    dssp = str_dssp=='True'
    sort_geo = str_sort_geo == 'True'
    geotag = geoInc.replace(':','_')
    if geotag != '':
        geotag = '_' + geotag
    csv_finalA = "Csv/PW_" + tagA + "_01_Geometry.csv"
    csv_finalB = "Csv/PW_" + tagB + "_01_Geometry.csv"
    if dssp:
        csv_finalA += "_dssp.csv"
        csv_finalB += "_dssp.csv"

    csv_correlations = "Csv/10_DivCorr_" + tagA + ".csv"
    html_filename = 'Html/20_ProteinCompare_' + tagA + '_' + tagB +geotag+str_sort_geo+ '.html'
    log_file = "Log/20_ProteinTop_dssp_" + tagA + tagB + geotag+ ".log"
    title = 'Williams Divergence from Trivial: Real versus Geometric: ' + tagA + ' vs ' + tagB
    if geoInc != '':
        title += ' Filtered on ' + geoInc
    openLog(log_file, csv_finalA)
    ###############################################################################################
    data = pd.read_csv(csv_finalA)
    dataB = pd.read_csv(csv_finalB)
    complete = pd.read_csv(csv_correlations)
    if sort_geo:
        complete = complete.sort_values(by='stat_ideal', ascending=False)
    else:
        complete = complete.sort_values(by='stat', ascending=False)
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
        data = data.query('occupancy == 1')
        data = data.query('bfactor <= 10')
        geos_to_abs = ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1']
        for gabs in geos_to_abs:
            data.loc[:,gabs] = abs(data.loc[:,gabs])
        aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
        data = data[data['aa'].isin(aa_list)]
        data = data[data['aa-1'].isin(aa_list)]
        data = data[data['aa+1'].isin(aa_list)]
        #remove omega 0  or pre omega 0
        data = data.query('`CA:C:N+1:CA+1` >= 100')
        data = data.query('`CA-1:C-1:N:CA` >= 100')
        if not dssp:
            data['ss'] = 'x'
            dataB['ss'] = 'x'

        for gabs in geos_to_abs:
            dataB[gabs] = abs(dataB[gabs])
        dataB = dataB.query('`CA:C:N+1:CA+1` >= 100')
        dataB = dataB.query('`CA-1:C-1:N:CA` >= 100')



    log(log_file, 'Create Williams Coefficient Maker')
    wcc = wcm.WilliamsDivergenceMaker(data, geos, density=density, log=0, norm=False, pval_iters=iters, delay_load=True)

    recreate_html = True
    if recreate_html:
        used_geos = []
        log(log_file, '### Creating html reports')
        rep_mak = hrm.HtmlReportMaker(title, html_filename, cols=6)
        rep_mak.addLineComment('Most correlated geos')
        i = -1
        while len(used_geos) < num_top:
            i += 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            statOrig = complete['stat'].values[i]
            statIdeal = complete['stat_ideal'].values[i]
            geoTagInc = True
            if geoInc != '':
                geoTagInc = False
                if geoA == geoInc or geoB == geoInc:
                    geoTagInc = True
            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos and geoTagInc:
                log(log_file, str(len(used_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                used_geos.append(geoA + '_' + geoB)
                statx, A, D, B = wcc.compareToObserved(dataB, geoA, geoB)
                maxV = max(np.max(A), np.max(B))
                rep_mak.addLineComment(geoA + ' | ' + geoB + ', Divergence from random=' + str(round(statOrig,4)) + ', Divergence from synthetic=' + str(round(statIdeal,4)))
                rep_mak.addPlot2d(data, 'scatter', title=str(round(statOrig, 3)) + ' orig', geo_x=geoA, geo_y=geoB, hue=geoA, palette='mako')
                rep_mak.addPlot2d(dataB, 'scatter', title='synthetic', geo_x=geoA, geo_y=geoB, hue=geoA, palette='mako')
                rep_mak.addSurface(A, 'Original Data', cmin=0, cmax=maxV, palette='Blues', colourbar=False)
                rep_mak.addSurface(D, 'Difference Data stat=' + str(round(statIdeal, 3)), cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                rep_mak.addSurface(B, 'Geometric Data', cmin=0, cmax=maxV, palette='Reds', colourbar=False)

        log(log_file, '############# Least correlated ###########################')
        rep_mak.addLineComment('Least correlated geos')
        i = len(complete.index)
        last_geos = []
        while len(last_geos) < num_bottom:
            i -= 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            statOrig = complete['stat'].values[i]
            statIdeal = complete['stat_ideal'].values[i]
            geoTagInc = True
            if geoInc != '':
                geoTagInc = False
                if geoA == geoInc or geoB == geoInc:
                    geoTagInc = True
            if geoA + '_' + geoB not in last_geos and geoB + '_' + geoA not in last_geos and geoTagInc:
                log(log_file, str(len(last_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                last_geos.append(geoA + '_' + geoB)
                stat, A, D, B = wcc.compareToObserved(dataB, geoA, geoB)
                maxV = max(np.max(A), np.max(B))
                rep_mak.addLineComment(geoA + ' | ' + geoB + ', Divergence from random=' + str(round(statOrig,4)) + ', Divergence from synthetic=' + str(round(statIdeal,4)))
                rep_mak.addPlot2d(data, 'scatter', title=str(round(statOrig, 3)) + ' orig', geo_x=geoA, geo_y=geoB,  hue=geoA, palette='mako')
                rep_mak.addPlot2d(dataB, 'scatter', title='synthetic', geo_x=geoA, geo_y=geoB, hue=geoA, palette='mako')
                rep_mak.addSurface(A, 'Original Data', cmin=0, cmax=maxV, palette='Blues', colourbar=False)
                rep_mak.addSurface(D, 'Difference Data stat=' + str(round(statIdeal, 3)), cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
                rep_mak.addSurface(B, 'Geometric Data', cmin=0, cmax=maxV, palette='Reds', colourbar=False)

        log(log_file, 'Finally print out to ' + html_filename)
        rep_mak.printReport()


###################################################################################
if __name__ == '__main__':
    globals()['proteinTopCompare'](sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])




