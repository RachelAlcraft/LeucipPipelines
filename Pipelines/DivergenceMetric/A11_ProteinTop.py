'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
modify_csv,recreate_divergence,recreate_html = True,True,True
csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
html_filename = 'Html/1a_ProteinTop.html'
title = 'Williams Divergence from Trivial: Most and Least Correlated'
num_top = 40
num_bottom = 10
iters = 1000
density = 1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help

from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsDivergenceMaker as wcm

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
###############################################################################################
data = pd.read_csv(csv_final)
geos = []
for col in data.columns:
    if ':' in col and 'info' not in col and '(' not in col and '{' not in col:
    #if ':' in col and 'info' not in col:
        geos.append(col)

#geos = ['N:O','N:CA','CA:C','N:CA:C:N+1']

if modify_csv:
    print('### Applying filters to csv data')
    geos_to_abs = ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1']
    for gabs in geos_to_abs:
        data[gabs] = abs(data[gabs])
    data['aa-1'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['infoN:N-1'], 1), axis=1)
    data['aa+1'] = data.apply(lambda row: help.applyDetailFromOther('AA', row['infoN:N+1'], 1), axis=1)

print('Create Williams Coefficient Maker')
wcc = wcm.WilliamsDivergenceMaker(data,geos,density=density,log=1,norm=False,pval_iters=iters,delay_load=(not recreate_divergence))

if recreate_divergence:
    complete = wcc.getCoefficientsDataFrame()
    complete = complete.sort_values(by='stat', ascending=False)
    complete.to_csv('Csv/11_Correlations.csv',index=False)
    print(complete)

complete = pd.read_csv('../WilliamsDivergence/Csv/11_Correlations.csv')
if recreate_html:
    used_geos = []
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title, html_filename, cols=6)
    rep_mak.addLineComment('Most correlated geos')
    i = -1
    while len(used_geos) < num_top:
        i += 1
        geoA = complete['geoA'].values[i]
        geoB = complete['geoB'].values[i]
        if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos:
            print(len(used_geos),geoA,geoB,'.........')
            used_geos.append(geoA + '_' + geoB)
            div = wcc.getCorrelation([geoA,geoB])
            cm_data = wcc.data
            df_rand = wcc.randomiseData(cm_data, [geoA, geoB])
            stat, pvalue, A, D, B = div.stat, div.p_value, div.histAB, div.diffAB, div.convAB
            mean, sd, hist = div.p_mean, div.p_std, div.p_hist
            maxV = max(np.max(A), np.max(B))
            rep_mak.addPlot2d(data, 'seaborn', title=str(round(stat, 3)) + ' orig', geo_x=geoA, geo_y=geoB, hue='aa+1',palette='Paired')
            rep_mak.addPlot2d(df_rand, 'scatter', title='rand', geo_x=geoA, geo_y=geoB, hue=geoA)
            if len(hist['divergence']) > 0:
                crit_val = round(wcc.getCriticalValue(geoA, geoB, 0.95), 3)
                rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence',  title='mean=' + str(round(mean, 3)) + ' sd=' + str(round(sd, 3)) + ' crit5%=' + str(crit_val),bins=50)
            else:
                rep_mak.addBoxComment(('No histogram calculated'))
            rep_mak.addSurface(A, 'Original Data', cmin=0, cmax=maxV, palette='Blues', colourbar=False)
            rep_mak.addSurface(D, 'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)),  cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
            rep_mak.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds', colourbar=False)

    print('############# Least correlated ###########################')
    rep_mak.addLineComment('Least correlated geos')
    i = len(complete.index)
    last_geos = []
    while len(last_geos) < num_bottom:
        i -= 1
        geoA = complete['geoA'].values[i]
        geoB = complete['geoB'].values[i]

        if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos:
            print(len(last_geos), geoA, geoB, '.........')
            last_geos.append(geoA + '_' + geoB)
            div = wcc.getCorrelation([geoA, geoB])
            cm_data = wcc.data
            df_rand = wcc.randomiseData(cm_data, [geoA, geoB])
            stat, pvalue, A, D, B = div.stat, div.p_value, div.histAB, div.diffAB, div.convAB
            mean, sd, hist = div.p_mean, div.p_std, div.p_hist
            maxV = max(np.max(A), np.max(B))
            rep_mak.addPlot2d(cm_data, 'seaborn', title=str(round(stat, 3)) + ' orig', geo_x=geoA, geo_y=geoB, hue='aa+1',palette='Paired')
            rep_mak.addPlot2d(df_rand, 'scatter', title='rand', geo_x=geoA, geo_y=geoB, hue=geoA)
            if len(hist['divergence']) > 0:
                crit_val = round(wcc.getCriticalValue(geoA, geoB, 0.95), 3)
                rep_mak.addPlot1d(hist, 'histogram', geo_x='divergence', title='mean=' + str(round(mean, 3)) + ' sd=' + str(round(sd, 3)) + ' crit5%=' + str(crit_val),bins=50)
            else:
                rep_mak.addBoxComment(('No histogram calculated'))
            rep_mak.addSurface(A, 'Original Data', cmin=0, cmax=maxV, palette='Blues', colourbar=False)
            rep_mak.addSurface(D, 'Difference Data stat=' + str(round(stat, 3)) + ' pvalue=' + str(round(pvalue, 3)), cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
            rep_mak.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds', colourbar=False)

    rep_mak.printReport()




