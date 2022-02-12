'''
Rachel Alcraft: 31/01/2022
Script using LeucipPy Protein Geometry Library
This script:
 1) Loads BioPython structures from pdb files in a directory
 2) Creates a dataframe that correlates chosen geoemtric measures for all those structures, you can use this to process the data elsewhere
 3) Displays chosen correlations of those geos in an html report driven by matplotlib and seaborn

 The leucipPy library can be installed from TestPyPi with the command

 pip install -i https://test.pypi.org/simple/ LeucipPy-pkg-RachelAlcraft

 There is a google colab which demonstrates it here

 https://colab.research.google.com/drive/1Ut5HXQQgE3sNuAMo7IVR3aqwUVOYDy12#scrollTo=2ubfgnaNP_Gs
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
import numpy as np

recreate_csv,modify_csv,recreate_html = False,True,True
# directory and file names
#directory of the pdb files
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data_redo/'
#full path, no path saves to current directory
csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
html_filename = 'Html/WilliamsCoeff.html'
#the geometric measures for geometry calculations

title = 'Williams Coefficient' # used at the header of the html report

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import DataFrameMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsCoefficientMaker as wcm
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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

#geos = ['N:CA','C:O','C:N+1','CA:C']
if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=2)
    print('Create Williams Coefficient Maker')
    wcc = wcm.WilliamsCoefficientMaker(data,geos,bins=20,log=1,norm=True)

    complete = wcc.getCoefficientsDataFrame()
    complete.to_csv('Csv/Correlations.csv')
    summary = complete.groupby(by='geoA').sum()
    summary = summary.sort_values(by='stat')
    summary.to_csv('Csv/Summary.csv')


    #fig, ax = plt.subplots()
    #splt = sns.barplot(x=summary['stat'].values,y=summary.index)
    #plt.xlabel('sum of coefficients')
    #plt.ylabel('geo')
    #rep_mak.addPlotOnly(fig, ax)

    # remove the low correlations
    rem_df = complete[['geoA','geoB','stat']]
    cut_off = 50
    for g in range(0,cut_off):
        geo = summary.index[g]
        rem_df = rem_df[rem_df['geoA'] != geo]
        rem_df = rem_df[rem_df['geoB'] != geo]
        #print(geo,rem_df.shape)

    #print(rem_df)

    red_piv = rem_df.pivot('geoA','geoB','stat')
    for i in range(0,len(red_piv.columns)):
        red_piv = red_piv.sort_values(by=red_piv.columns[i], ascending=False)
        red_piv = red_piv[list(red_piv.index.values)]
    for i in range(len(red_piv.columns)-1,-1,-1):
        red_piv = red_piv.sort_values(by=red_piv.columns[i], ascending=False)
        red_piv = red_piv[list(red_piv.index.values)]

    fig, ax = plt.subplots()
    sns.heatmap(red_piv, annot=False, fmt='.2f', linewidth=.5, ax=ax, cmap='inferno_r', mask=red_piv.isnull(), vmin=0, xticklabels=red_piv.columns,yticklabels=red_piv.index)
    plt.title('Heatmap with least correlated removed')
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6, rotation=90)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=6, rotation=0)
    ax.tick_params(which="both", bottom=True)
    rep_mak.addPlotOnly(fig, ax)

    #for frac,asc,ttl in [[1,False,'Heatmap of all geos'],[0.2,False,'Heatmap top 20%'],[0.2,True,'Heatmap bottom 20%']]:
    heat=wcc.getCoefficientsDataFrame(as_pivot=True,fraction=1,asc=False)

    for i in range(0, len(heat.columns)):
        heat = heat.sort_values(by=heat.columns[i], ascending=False)
        heat = heat[list(heat.index.values)]
    for i in range(len(heat.columns) - 1, -1, -1):
        heat = heat.sort_values(by=heat.columns[i], ascending=False)
        heat = heat[list(heat.index.values)]


    fig, ax = plt.subplots()
    sns.heatmap(heat, annot=False, fmt='.2f', linewidth=.5, ax=ax, cmap='inferno_r',mask=heat.isnull(),vmin=0)
    plt.title('Heatmap of all geos')
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6,rotation=90)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=6, rotation=0)
    rep_mak.addPlotOnly(fig, ax)


    print('Make scatters')
    complete = complete.sort_values(by='stat',ascending=False)

    rep_mak.changeColNumber(4)
    max_plots = 50
    for i in range(0,max_plots):
    #for i in range(0,len((summary.index))):
        geoA = complete['geoA'].values[i]
        geoB = complete['geoB'].values[i]
        stat = complete['stat'].values[i]
        print('...', geoA,geoB,stat)
        #stat = summary['stat'].values[i]
        #heat = wcc.getRelativePlotCoefficientsDataFrame(geo,as_pivot=True)
        #fig, ax = plt.subplots()
        #plot = sns.heatmap(heat, annot=False, fmt='.2f', linewidth=.5, ax=ax, cmap='inferno_r', mask=heat.isnull(),vmin=0,vmax=1)
        #ttl = geo + ' ' + str(round(stat, 2))
        #plt.title(ttl)
        #ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6, rotation=90)
        #ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=6, rotation=0)
        #print(geo, stat)
        #rep_mak.addPlotOnly(fig, ax)
        #for j in range(i+1, max_plots):
        #    geoB = summary.index[j]

        rep_mak.addPlot2d(data, 'scatter', title=geoA + '|' + geoB, geo_x=geoA, geo_y=geoB, hue=geoA)

        print('Scatter', geoA, geoB)
        geo, stat, A, D, B = wcc.getCorrelation([geoA, geoB])
        maxV = max(np.max(A), np.max(B))
        rep_mak.addSurface(A, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Blues', colourbar=False)
        rep_mak.addSurface(D, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
        rep_mak.addSurface(B, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Reds', colourbar=False)

    rep_mak.addLineComment(('Least correlated 10'))
    for i in range(len(complete.index)-10,len(complete.index)):
        geoA = complete['geoA'].values[i]
        geoB = complete['geoB'].values[i]
        stat = complete['stat'].values[i]
        print('...', geoA,geoB,stat)

        rep_mak.addPlot2d(data, 'scatter', title=geoA + '|' + geoB, geo_x=geoA, geo_y=geoB, hue=geoA)

        print('Scatter', geoA, geoB)
        geo, stat, A, D, B = wcc.getCorrelation([geoA, geoB])
        maxV = max(np.max(A), np.max(B))
        rep_mak.addSurface(A, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Blues', colourbar=False)
        rep_mak.addSurface(D, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=-1 * maxV, cmax=maxV, palette='RdBu', colourbar=False)
        rep_mak.addSurface(B, geoA + '|' + geoB + ' bins =' + str(wcc.bins) + ' coeff=' + str(round(stat, 3)),cmin=0 * maxV, cmax=maxV, palette='Reds', colourbar=False)



    rep_mak.printReport()




