'''
Rachel Alcraft: 31/01/2022
Script using LeucipPy Protein Geometry Library
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
reference_geo = 'CA-1:CA:CA+1'
tag = 'CA3'
recreate_csv,modify_csv,recreate_html = False,True,True
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data_redo/'
csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
html_filename = 'Html/WilliamsCoeff_GLY_' + tag + '.html'
title = 'Williams Coefficient GLY ' + tag# used at the header of the html report

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsCoefficientMaker as wcm
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
###############################################################################################
data = pd.read_csv(csv_final)
geos = []
for col in data.columns:
    if ':' in col and 'info' not in col and '(' not in col and '{' not in col:
    #if ':' in col and 'info' not in col:
        geos.append(col)
if modify_csv:
    print('### Applying filters to csv data')
    geos_to_abs = ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1']
    for gabs in geos_to_abs:
        data[gabs] = abs(data[gabs])

if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=2)
    print('Create Williams Coefficient Maker')
    wcc = wcm.WilliamsCoefficientMaker(data,geos,bins=20,log=1,norm=True,cut_off=25)

    complete = wcc.getCoefficientsDataFrame()
    complete.to_csv('Csv/Correlations.csv')
    summary = complete.groupby(by='geoA').sum()
    summary = summary.sort_values(by='stat')
    summary.to_csv('Csv/Summary.csv')

    # remove all but geoA == TAU
    rem_df = complete[['geoA','geoB','stat']]
    rem_df = rem_df[rem_df['geoA'] == reference_geo]

    print('Make scatters')
    rem_df = rem_df.sort_values(by='stat',ascending=False)

    rep_mak.changeColNumber(4)
    max_plots = 50
    for i in range(0,max_plots):
        geoA = rem_df['geoA'].values[i]
        geoB = rem_df['geoB'].values[i]
        stat = rem_df['stat'].values[i]
        print('...', geoA,geoB,stat)
        title = geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3))
        rep_mak.addLineComment(title)
        rep_mak.addPlot2d(data, 'scatter', title=title, geo_x=geoA, geo_y=geoB, hue='N:CA:C',palette='Spectral')
        print('Scatter', geoA, geoB)
        geo, stat, A, D, B = wcc.getCorrelation([geoA, geoB])
        maxV = max(np.max(A), np.max(B))
        rep_mak.addSurface(A, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Blues', colourbar=False)
        rep_mak.addSurface(D, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
        rep_mak.addSurface(B, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Reds', colourbar=False)

    rep_mak.addLineComment(('Least correlated 10'))
    for i in range(len(rem_df.index)-10,len(rem_df.index)):
        geoA = rem_df['geoA'].values[i]
        geoB = rem_df['geoB'].values[i]
        stat = rem_df['stat'].values[i]
        print('...', geoA,geoB,stat)
        title = geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3))
        rep_mak.addLineComment(title)
        rep_mak.addPlot2d(data, 'scatter', title=title, geo_x=geoA, geo_y=geoB, hue=geoA)

        print('Scatter', geoA, geoB)
        geo, stat, A, D, B = wcc.getCorrelation([geoA, geoB])
        maxV = max(np.max(A), np.max(B))
        rep_mak.addSurface(A, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Blues', colourbar=False)
        rep_mak.addSurface(D, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
        rep_mak.addSurface(B, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Reds', colourbar=False)



    rep_mak.printReport()




