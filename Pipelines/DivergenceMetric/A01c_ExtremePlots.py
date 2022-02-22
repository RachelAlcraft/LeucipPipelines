

import math
import random
import sys

from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsDivergenceMaker as wcm
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from datetime import datetime as dt
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

runs = []
runs.append(['spot',500])
runs.append(['vert',500])
runs.append(['line',500])
runs.append(['covar',500])
runs.append(['rand',500])

html_file = "Html/A01c_Baseline_Extremes.html"
csv_file = "Csv/A01c_Baseline_Extremes.csv"
log_file = "Log/A01c_Baseline_Extremes.log"
rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Baseline",html_file,cols=2)
openLog(log_file,'')

for randOrline,iters in runs:
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    fake_geos =['geoA','geoB']
    samplesize = 1000
    bins = 5
    density = samplesize/(bins*bins)
    normed = True

    log(log_file,'Create samples')
    dic_sample_bin = {}
    dic_fake_1 = {'geoA': [], 'geoB': []}
    for i in range(0,samplesize):
        if randOrline == 'spot':
            dic_fake_1['geoA'].append(50)
            dic_fake_1['geoB'].append(50)
        elif randOrline == 'line':
            dic_fake_1['geoA'].append(i)
            dic_fake_1['geoB'].append(i)
        elif randOrline == 'vert':
            dic_fake_1['geoA'].append(50)
            dic_fake_1['geoB'].append(i)
        elif randOrline == 'covar':
            l = i%20
            vi = 0.3
            count = 10#samplesize
            dic_fake_1['geoA'].append(np.random.normal(l, int(count * vi)))
            dic_fake_1['geoB'].append(np.random.normal(l, int(count * vi)))
        elif randOrline == 'rand':
            l = 20
            count = 10#samplesize
            dic_fake_1['geoA'].append(random.randint(0, l))
            dic_fake_1['geoB'].append(random.randint(0, l))

    df_sample = pd.DataFrame.from_dict(dic_fake_1)
    print(df_sample)

    log(log_file, 'samples=' + str(samplesize) + ' bins=' + str(density))
    cm =wcm.WilliamsDivergenceMaker(df_sample,fake_geos,density=density,log=2,norm=normed,pval_iters=iters,delay_load=False)

    #add a seaborn lineplot where there is NO sin or cos involved
    log(log_file, 'Make distribution plots')
    rep.addLineComment('Plots of distributions')

    log(log_file, 'Plots bins=' + str(density) + ' size=' + str(samplesize))
    cm_data = cm.data
    print(cm_data)
    df_rand = cm.randomiseData(cm_data,['geoA', 'geoB'])
    div = cm.getCorrelation(['geoA', 'geoB'])
    print(div.convAB)
    stat,pvalue,A,D,B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
    mean,sd,hist = div.p_mean,div.p_std,div.p_hist
    maxV = max(np.max(A),np.max(B))
    rep.changeColNumber(8)
    rep.addPlot2d(df_sample, 'scatter', title='Observed Data stat=' + str(round(stat,3)) + ' pvalue=' + str(round(pvalue,3)), geo_x='geoA', geo_y='geoB', hue='geoA')
    rep.addPlot2d(df_rand, 'scatter',title='rand, size=' + str(samplesize), geo_x='geoA', geo_y='geoB', hue='geoA')
    rep.addPlot1d(df_sample,'histogram','geoA')
    rep.addPlot1d(df_sample, 'histogram', 'geoB')
    rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues',colourbar=False)
    rep.addSurface(D, 'Difference Data stat=' + str(round(stat,3)) + ' pvalue=' + str(round(pvalue,3)), cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
    rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds',colourbar=False)
    if len(hist['divergence']) > 0:
        crit_val = round(cm.getCriticalValue('geoA', 'geoB', 0.95), 3)
        rep.addPlot1d(hist, 'histogram', geo_x='divergence', title='mean=' + str(round(mean, 3)) + ' sd=' + str(round(sd, 3)) + ' crit5%=' + str(crit_val),bins=50)
    else:
        rep.addBoxComment(('No histogram calculated'))

########################
log(log_file, 'Finally print out')
rep.printReport()


