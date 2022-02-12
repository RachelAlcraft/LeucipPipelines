

import math
import random
import sys

import scipy
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import DataFrameMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsCoefficientMaker as wcm
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def syntheticRunner(normed_corrstr,sample_sizestr):
    normed_corr = normed_corrstr=='True'
    sample_size = int(sample_sizestr)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    rep = hrm.HtmlReportMaker("Williams Divergence From Trivial Coefficient","Html/Synthetic_"+str(normed_corr) + str(sample_size) +".html",cols=8)
    dic_fake = {'Line':[],'Const':[],'Rand':[],'Sina':[],'Cosa':[],'Norm':[],'Sinb':[],'Cosb':[]}
    fake_geos =['Line','Const','Rand','Norm','Sina','Cosa','Sinb','Cosb']
    long_freq = sample_size/100
    short_freq = sample_size/10
    for i in range(0,sample_size):
        dic_fake['Line'].append(i)
        dic_fake['Const'].append(5)
        dic_fake['Rand'].append(random.randint(0,5000))
        dic_fake['Sina'].append(math.sin(i/long_freq))
        dic_fake['Cosa'].append(math.cos(i/long_freq))
        dic_fake['Norm'].append(random.randint(0,5000))
        #dic_fake['Norm'].append(np.random.normal(0, 1))
        dic_fake['Sinb'].append(math.sin(i/short_freq))
        dic_fake['Cosb'].append(math.cos(i/short_freq))

    df_fake = pd.DataFrame.from_dict(dic_fake)
    bins_all = [5,10,20,25,40,50,60,80,100,150,200]
    bins_red = [5,10,20,50,100]
    dic_cm = {}
    for bin in bins_all:
        dic_cm[bin] =wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=bin,log=1,norm=normed_corr)

    print('Creating histograms')
    for bin in [20]:
        cm = dic_cm[bin]
        rep.addLineComment('Williams Coefficient Maker, bins=' + str(cm.bins) +  ' samples=' + str(sample_size))
        for geo in fake_geos:
            stat = cm.correlations1d[geo]
            rep.addPlot1d(df_fake,'histogram',geo,hue='',bins=50,title=str(round(stat,3)))


    rep.changeColNumber(4)
    print('Compare coefficients')
    #rep.addBoxComment('Comparing the change in coefficients')
    lines_dic_all_sina = {'bins':[],'stat':[],'corr':[],'rank':[],'change':[]}
    lines_dic_all_sinb = {'bins':[],'stat':[],'corr':[],'rank':[],'change':[]}
    lines_dic_all_line = {'bins':[],'stat':[],'corr':[],'rank':[],'change':[]}
    lines_dic_all_none = {'bins':[],'stat':[],'corr':[],'rank':[],'change':[]}
    for i in range(0,len(fake_geos)):
        geoA = fake_geos[i]
        for j in range(i+1,len(fake_geos)):
            geoB = fake_geos[j]
            if geoA != geoB:
                last_stat = 0
                for bin in bins_all:
                    cm = dic_cm[bin]
                    print('Coeff',geoA,geoB)
                    geo,stat,A,D,B = cm.getCorrelation([geoA, geoB])
                    #stat = stat / (sample_size/bin**2)
                    change = 0
                    if last_stat != 0:
                        change = last_stat - stat
                    if stat > 0.001:
                        tau, p = scipy.stats.kendalltau(df_fake[geoA], df_fake[geoB])
                    #make dataframes
                        #if 'Sin' in geoA + geoB and 'Cos' in geoA+geoB or (geoA+geoB == 'SinaSinb' or geoA+geoB == 'CosaCosb') or 'LineSin' in geoA+geoB or 'LineCos' in geoA+geoB:
                        if ('Sina' in geoA+geoB or 'Cosa' in geoA + geoB):# and ('Rand' not in geoA+geoB) and ('Norm' not in geoA+geoB):
                            lines_dic_all_sina['corr'].append(geoA + geoB)
                            lines_dic_all_sina['bins'].append(cm.bins)
                            lines_dic_all_sina['stat'].append(stat)
                            lines_dic_all_sina['rank'].append(tau)
                            lines_dic_all_sina['change'].append(change)
                        if 'Sinb' in geoA+geoB or 'Cosb' in geoA + geoB:# and ('Rand' not in geoA+geoB) and ('Norm' not in geoA+geoB):
                            lines_dic_all_sinb['corr'].append(geoA + geoB)
                            lines_dic_all_sinb['bins'].append(cm.bins)
                            lines_dic_all_sinb['stat'].append(stat)
                            lines_dic_all_sinb['rank'].append(tau)
                            lines_dic_all_sinb['change'].append(change)
                        if 'Line' in geoA+geoB:# or 'Rand' in geoA + geoB or 'Norm' in geoA + geoB:
                            lines_dic_all_line['corr'].append(geoA + geoB)
                            lines_dic_all_line['bins'].append(cm.bins)
                            lines_dic_all_line['stat'].append(stat)
                            lines_dic_all_line['rank'].append(tau)
                            lines_dic_all_line['change'].append(change)
                        if 'Rand' in geoA + geoB or 'Norm' in geoA + geoB:
                            lines_dic_all_none['corr'].append(geoA + geoB)
                            lines_dic_all_none['bins'].append(cm.bins)
                            lines_dic_all_none['stat'].append(stat)
                            lines_dic_all_none['rank'].append(tau)
                            lines_dic_all_none['change'].append(change)
                        last_stat = stat

    lines_df_sina = pd.DataFrame.from_dict(lines_dic_all_sina).round(3).sort_values(by='stat',ascending=False)
    lines_df_sinb = pd.DataFrame.from_dict(lines_dic_all_sinb).round(3).sort_values(by='stat',ascending=False)
    lines_df_line = pd.DataFrame.from_dict(lines_dic_all_line).round(3).sort_values(by='stat',ascending=False)
    lines_df_none = pd.DataFrame.from_dict(lines_dic_all_none).round(3).sort_values(by='stat',ascending=False)
    print(lines_df_none)

    #add a seaborn lineplot where there is a sin or cos involved
    rep.addLineComment('Compare coefficients')
    y_val = 'stat'
    ylim = 1
    if normed_corr:
        ylim=10
    fig, ax = plt.subplots()
    sns.lineplot(data=lines_dic_all_sina,x='bins',y=y_val,hue='corr')#,palette='Paired')
    plt.title('Change in correlations over bins, Sina+Cosa')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,ylim)
    rep.addPlotOnly(fig, ax)
    #add a seaborn lineplot where there is NO sin or cos involved
    fig, ax = plt.subplots()
    sns.lineplot(data=lines_dic_all_sinb,x='bins',y=y_val,hue='corr')#,palette='Paired')
    plt.title('Change in correlations over bins, Sinb+Cosb')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,ylim)
    rep.addPlotOnly(fig, ax)
    #add a seaborn lineplot where there is NO sin or cos involved
    fig, ax = plt.subplots()
    sns.lineplot(data=lines_dic_all_line,x='bins',y=y_val,hue='corr')#,palette='Paired')
    plt.title('Change in correlations over bins, Line')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,ylim)
    rep.addPlotOnly(fig, ax)
    #add a seaborn lineplot where there is NO sin or cos involved
    fig, ax = plt.subplots()
    sns.lineplot(data=lines_dic_all_none,x='bins',y=y_val,hue='corr')#,palette='Paired')
    plt.title('Change in correlations over bins, Norm+Rand')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,ylim)
    rep.addPlotOnly(fig, ax)


    rep.addLineComment('Data in the dataframes, bins=' + str(10))
    rep.addDataFrame(lines_df_sina[lines_df_sina['bins']==10],title='Short sine waves')
    rep.addDataFrame(lines_df_sinb[lines_df_sinb['bins']==10],title='Long sine waves')
    rep.addDataFrame(lines_df_line[lines_df_line['bins']==10],title='Correlated with a line')
    rep.addDataFrame(lines_df_none[lines_df_none['bins']==10],title='Random distributions')
    rep.addLineComment('Data in the dataframes, bins=' + str(25))
    rep.addDataFrame(lines_df_sina[lines_df_sina['bins']==25],title='Short sine waves')
    rep.addDataFrame(lines_df_sinb[lines_df_sinb['bins']==25],title='Long sine waves')
    rep.addDataFrame(lines_df_line[lines_df_line['bins']==25],title='Correlated with a line')
    rep.addDataFrame(lines_df_none[lines_df_none['bins']==25],title='Random distributions')



    rep.addLineComment('Plots of changing correlations')
    rep.changeColNumber(6)
    print('Show some coefficients')
    lines_dic = {}
    lines_dic['bins'] = [5,10,20,50,100]
    for i in range(0,len(fake_geos)):
        geoA = fake_geos[i]
        for j in range(i+1,len(fake_geos)):
            geoB = fake_geos[j]
            if geoA != geoB:
                lines_dic[geoA + geoB] = []
                rep.addPlot2d(df_fake, 'scatter',title=geoA + '|' + geoB , geo_x=geoA, geo_y=geoB, hue=geoA)
                for bin in bins_red:
                    cm = dic_cm[bin]
                    print('Scatter',geoA,geoB)
                    geo,stat,A,D,B = cm.getCorrelation([geoA, geoB])
                    maxV = max(np.max(A),np.max(B))
                    rep.addSurface(D, geoA + '|' + geoB +  ' bins =' + str(cm.bins) + ' coeff='+str(round(stat,3)), cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
                    lines_dic[geoA + geoB].append(stat)

                lines_df = pd.DataFrame.from_dict(lines_dic)
                #rep.addPlot2d(lines_df, 'scatter', title=geoA + '|' + geoB, geo_x='bins', geo_y=geoA + geoB, hue='bins',yrange=[0,1])




    rep.addLineComment('Scatters')
    print('Creating scatters')
    rep.changeColNumber(5)
    for i in range(0,len(fake_geos)):
        geoA = fake_geos[i]
        for j in range(i+1,len(fake_geos)):
            geoB = fake_geos[j]
            if geoA != geoB:

                print('Scatter',geoA,geoB)
                cm = dic_cm[20]
                df_rand = cm.randomiseData(df_fake,[geoA,geoB])
                geo,stat,A,D,B = cm.getCorrelation([geoA, geoB])
                if stat > 0.001:
                    maxV = max(np.max(A),np.max(B))
                    rep.addLineComment(geoA + ' ' + geoB + ' wdtc=' + str(round(stat,3)))
                    rep.addPlot2d(df_fake, 'scatter', title=str(round(stat,3)),geo_x=geoA, geo_y=geoB, hue=geoA)
                    rep.addPlot2d(df_rand, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA)
                    rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues',colourbar=False)
                    rep.addSurface(D, 'Difference Data', cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
                    rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds',colourbar=False)


    rep.printReport()



if __name__ == '__main__':
    globals()['syntheticRunner'](sys.argv[1],sys.argv[2])