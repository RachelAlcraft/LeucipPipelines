

import math
import random

from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsDivergenceMaker as wcm
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
normed_corr = False
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

runs = []
#runs.append(['rand',5])
runs.append(['rand',500])
runs.append(['line',500])
runs.append(['covar',500])

for randOrline,sample_size in runs:
    print('Starting loop',randOrline,sample_size)
    html_file = "Html/02_Normalised_" + randOrline + str(sample_size) + ".html"
    print(html_file)
    csv_file = "Csv/01_Baseline_" + randOrline + str(sample_size) + ".csv"
    #csv_file = "Csv/01_stats" + randOrline + ".csv"

    rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Baseline",html_file,cols=2)
    fake_geos =['RandA','RandB']
    df_fakea = pd.read_csv(csv_file)
    dic_fake = {}
    dic_fake['p_value'] =df_fakea['p_value'].values
    dic_fake['bins'] =df_fakea['bins'].values
    dic_fake['samples'] =df_fakea['samples'].values
    dic_fake['stat'] =df_fakea['stat'].values
    dic_fake['stat2'] =[]
    dic_fake['density'] =[]
    for i in range(0,len(df_fakea['stat'].values)):
        stat = df_fakea['stat'].values[i]
        samp = df_fakea['samples'].values[i]
        bn = df_fakea['bins'].values[i]
        density = samp/(bn*bn)
        dic_fake['density'].append(density)
        #stat2 = (math.log(stat*density))
        #stat2 = ((stat * density))
        stat2 = (math.log(stat * samp/bn))
        dic_fake['stat2'].append(stat2)

    df_fake = pd.DataFrame.from_dict(dic_fake).round(4)

    lineplots = []
    lineplots.append(['bins','stat','samples','Change in correlations all bins',True,False,'',0])
    lineplots.append(['bins', 'stat', 'samples', '<=100 bins', True,True,'bins',100])
    lineplots.append(['bins', 'stat2', 'samples', 'stat2', False, False, 'bins', 100])
    lineplots.append(['bins', 'stat2', 'samples', 'stat2', False, True, 'bins', 100])
    lineplots.append(['bins', 'p_value', 'samples', '', False, False, 'bins', 100])
    lineplots.append(['bins', 'p_value', 'samples', '', False, True, 'bins', 100])
    lineplots.append(['density', 'stat2', 'samples', '', False, False, 'bins', 100])
    lineplots.append(['density', 'stat2', 'samples', '', False, True, 'density', 1000])
    lineplots.append(['density', 'stat2', 'bins', '', False, False, 'bins', 100])
    lineplots.append(['density', 'stat2', 'bins', '', False, True, 'density', 1000])
    lineplots.append(['samples', 'stat2', 'density', '', False, True, 'density', -100])
    lineplots.append(['samples', 'stat2', 'bins', '', False, False, 'bins', 100])


    for xx,yy,hh,tt,lim,change,col,val in lineplots:
        fig, ax = plt.subplots()
        df_use = df_fake
        if change:
            if val > 0:
                df_use = df_fake[df_fake[col] <= val]
            else:
                print(col,val)
                df_use = df_fake[df_fake[col] >= val]

            sns.lineplot(data=df_use, x=xx, y=yy, hue=hh, palette='viridis')
        else:
            sns.lineplot(data=df_use,x=xx,y=yy,hue=hh,palette='Paired')
        plt.title(tt)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        if lim:
            plt.ylim(0,1)
        rep.addPlotOnly(fig, ax)


    #add a seaborn lineplot where there is NO sin or cos involved
    rep.addLineComment('Data in the dataframes')
    for bin in df_fake['bins'].unique():
        rep.addDataFrame(df_fake[df_fake['bins']==bin],title='bin size=' + str(bin))

    rep.addPlot2d(df_fake,'scatter',geo_x='bins',geo_y='stat2',hue='samples')
    rep.addPlot2d(df_fake, 'scatter', geo_x='samples', geo_y='stat2', hue='bins')
    #rep.addPlot2d(df_fake, 'scatter', geo_x='density', geo_y='stat2', hue='bins')
    rep.addPlot2d(df_fake, 'scatter', geo_x='samples', geo_y='bins', hue='stat2')
    rep.addPlot2d(df_fake, 'seaborn', geo_x='samples', geo_y='bins', hue='density')



    rep.printReport()


