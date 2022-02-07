

import math
import random

from LeucipPy import BioPythonMaker as bpm
from LeucipPy import DataFrameMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsCoefficientMaker as wcm
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

rep = hrm.HtmlReportMaker("Williams Divergence From Trivial Coefficient","Html/synthetic.html",cols=8)
################ FALSE DATA #########################
normed_corr = True
dic_fake = {'Line':[],'Const':[],'Rand':[],'Sina':[],'Cosa':[],'Norm':[],'Sinb':[],'Cosb':[]}
for i in range(0,5000):
    dic_fake['Line'].append(i)
    dic_fake['Const'].append(5)
    dic_fake['Rand'].append(random.randint(0,5000))
    dic_fake['Sina'].append(math.sin(i/50))
    dic_fake['Cosa'].append(math.cos(i/50))
    dic_fake['Norm'].append(np.random.normal(0,1))
    dic_fake['Sinb'].append(math.sin(i/500))
    dic_fake['Cosb'].append(math.cos(i/500))


fake_geos =['Line','Const','Rand','Sina','Cosa','Norm','Sinb','Cosb']
df_fake = pd.DataFrame.from_dict(dic_fake)
cm_5 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=5,log=1,norm=normed_corr)
cm_10 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=10,log=1,norm=normed_corr)
cm_20 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=20,log=1,norm=normed_corr)
cm_25 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=25,log=1,norm=normed_corr)
cm_40 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=40,log=1,norm=normed_corr)
cm_50 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=50,log=1,norm=normed_corr)
cm_60 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=60,log=1,norm=normed_corr)
cm_80 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=80,log=1,norm=normed_corr)
cm_100 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=100,log=1,norm=normed_corr)
cm_150 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=150,log=1,norm=normed_corr)
cm_200 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=200,log=1,norm=normed_corr)
#cm_500 = wcm.WilliamsCoefficientMaker(df_fake,fake_geos,bins=500,log=1,norm=normed_corr)
print('Creating histograms')
cms = [cm_5,cm_10,cm_20,cm_50,cm_100]
cms_all = [cm_5,cm_10,cm_20,cm_25,cm_40,cm_50,cm_60,cm_80,cm_100,cm_150,cm_200]#,cm_500]
for cm in [cm_20]:
    rep.addLineComment('Williams Coefficient Maker, bins=' + str(cm.bins))
    for geo in fake_geos:
        stat = cm.correlations1d[geo]
        rep.addPlot1d(df_fake,'histogram',geo,hue='',bins=50,title=str(round(stat,3)))


rep.addLineComment('Compare coefficients')
rep.changeColNumber(2)
print('Compare coefficients')
#rep.addBoxComment('Comparing the change in coefficients')
lines_dic_all_sina = {'bins':[],'stat':[],'corr':[]}
lines_dic_all_sinb = {'bins':[],'stat':[],'corr':[]}
lines_dic_all_line = {'bins':[],'stat':[],'corr':[]}
lines_dic_all_none = {'bins':[],'stat':[],'corr':[]}
for i in range(0,len(fake_geos)):
    geoA = fake_geos[i]
    for j in range(i+1,len(fake_geos)):
        geoB = fake_geos[j]
        if geoA != geoB:
            for cm in cms_all:
                print('Coeff',geoA,geoB)
                geo,stat,A,D,B = cm.getCorrelation([geoA, geoB])
                if stat > 0.001:
                #make dataframes
                    #if 'Sin' in geoA + geoB and 'Cos' in geoA+geoB or (geoA+geoB == 'SinaSinb' or geoA+geoB == 'CosaCosb') or 'LineSin' in geoA+geoB or 'LineCos' in geoA+geoB:
                    if ('Sina' in geoA+geoB or 'Cosa' in geoA + geoB):# and ('Rand' not in geoA+geoB) and ('Norm' not in geoA+geoB):
                        lines_dic_all_sina['corr'].append(geoA + geoB)
                        lines_dic_all_sina['bins'].append(cm.bins)
                        lines_dic_all_sina['stat'].append(stat)
                    if 'Sinb' in geoA+geoB or 'Cosb' in geoA + geoB:# and ('Rand' not in geoA+geoB) and ('Norm' not in geoA+geoB):
                        lines_dic_all_sinb['corr'].append(geoA + geoB)
                        lines_dic_all_sinb['bins'].append(cm.bins)
                        lines_dic_all_sinb['stat'].append(stat)
                    if 'Line' in geoA+geoB:# or 'Rand' in geoA + geoB or 'Norm' in geoA + geoB:
                        lines_dic_all_line['corr'].append(geoA + geoB)
                        lines_dic_all_line['bins'].append(cm.bins)
                        lines_dic_all_line['stat'].append(stat)
                    if 'Rand' in geoA + geoB or 'Norm' in geoA + geoB:
                        lines_dic_all_none['corr'].append(geoA + geoB)
                        lines_dic_all_none['bins'].append(cm.bins)
                        lines_dic_all_none['stat'].append(stat)

lines_df_sina = pd.DataFrame.from_dict(lines_dic_all_sina)
lines_df_sinb = pd.DataFrame.from_dict(lines_dic_all_sinb)
lines_df_line = pd.DataFrame.from_dict(lines_dic_all_line)
lines_df_none = pd.DataFrame.from_dict(lines_dic_all_none)
#add a seaborn lineplot where there is a sin or cos involved
y_val = 'stat'
fig, ax = plt.subplots()
sns.lineplot(data=lines_dic_all_sina,x='bins',y=y_val,hue='corr')#,palette='Paired')
plt.title('Change in correlations over bins, Sina+Cosa')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylim(0,10)
rep.addPlotOnly(fig, ax)
#add a seaborn lineplot where there is NO sin or cos involved
fig, ax = plt.subplots()
sns.lineplot(data=lines_dic_all_sinb,x='bins',y=y_val,hue='corr')#,palette='Paired')
plt.title('Change in correlations over bins, Sinb+Cosb')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylim(0,10)
rep.addPlotOnly(fig, ax)
#add a seaborn lineplot where there is NO sin or cos involved
fig, ax = plt.subplots()
sns.lineplot(data=lines_dic_all_line,x='bins',y=y_val,hue='corr')#,palette='Paired')
plt.title('Change in correlations over bins, Line')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylim(0,10)
rep.addPlotOnly(fig, ax)
#add a seaborn lineplot where there is NO sin or cos involved
fig, ax = plt.subplots()
sns.lineplot(data=lines_dic_all_none,x='bins',y=y_val,hue='corr')#,palette='Paired')
plt.title('Change in correlations over bins, Norm+Rand')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylim(0,10)
rep.addPlotOnly(fig, ax)


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
            for cm in cms:
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
            df_rand = cm_20.randomiseData(df_fake,[geoA,geoB])
            geo,stat,A,D,B = cm_20.getCorrelation([geoA, geoB])
            if stat > 0.001:
                maxV = max(np.max(A),np.max(B))
                rep.addLineComment(geoA + ' ' + geoB + ' wdtc=' + str(round(stat,3)))
                rep.addPlot2d(df_fake, 'scatter', title=str(round(stat,3)),geo_x=geoA, geo_y=geoB, hue=geoA)
                rep.addPlot2d(df_rand, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA)
                rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues',colourbar=False)
                rep.addSurface(D, 'Difference Data', cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
                rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds',colourbar=False)


rep.printReport()



