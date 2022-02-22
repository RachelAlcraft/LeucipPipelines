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
    csv_correlations = "Csv/10_DivCorr_" + tag + ".csv"
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
        data = data.query('occupancy == 1')
        data = data.query('bfactor <= 10')
        geos_to_abs = ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1']
        for gabs in geos_to_abs:
            data[gabs] = abs(data[gabs])
        aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
        data = data[data['aa'].isin(aa_list)]
        data = data[data['aa-1'].isin(aa_list)]
        data = data[data['aa+1'].isin(aa_list)]


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
                df_rand = wcc.randomiseData(cm_data, [geoA, geoB])
                stat, pvalue, A, D, B = div.stat, div.p_value, div.histAB, div.diffAB, div.convAB
                mean, sd, hist = div.p_mean, div.p_std, div.p_hist
                maxV = max(np.max(A), np.max(B))
                rep_mak.addPlot2d(cm_data, 'seaborn', title='Observed Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), geo_x=geoA, geo_y=geoB, hue='aa',palette='tab20')
                rep_mak.addPlot2d(df_rand, 'scatter', title='rand', geo_x=geoA, geo_y=geoB, hue=geoA)
                if len(hist['divergence']) > 0:
                    crit_val = round(wcc.getCriticalValue(geoA, geoB, 0.95), 3)
                    rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence',  title='mean=' + str(round(mean, 3)) + ' sd=' + str(round(sd, 3)) + ' crit5%=' + str(crit_val),bins=50)
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
                df_rand = wcc.randomiseData(cm_data, [geoA, geoB])
                stat, pvalue, A, D, B = div.stat, div.p_value, div.histAB, div.diffAB, div.convAB
                mean, sd, hist = div.p_mean, div.p_std, div.p_hist
                maxV = max(np.max(A), np.max(B))
                rep_mak.addPlot2d(cm_data, 'seaborn', title='Observed Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), geo_x=geoA, geo_y=geoB, hue='aa',palette='tab20_r')
                rep_mak.addPlot2d(df_rand, 'scatter', title='rand', geo_x=geoA, geo_y=geoB, hue=geoA)
                if len(hist['divergence']) > 0:
                    crit_val = round(wcc.getCriticalValue(geoA, geoB, 0.95), 3)
                    rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence', title='mean=' + str(round(mean, 3)) + ' sd=' + str(round(sd, 3)) + ' crit5%=' + str(crit_val),bins=50)
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




