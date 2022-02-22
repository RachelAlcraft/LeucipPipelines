

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

def randomBaseline(randOrline,str_iters):
    normed_corr = False
    iters = int(str_iters)
    #randOrline = 'line' #rand or line or covar
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    html_file = "Html/X01_Baseline_"+randOrline+str_iters+".html"
    csv_file = "Csv/X01_Baseline_"+randOrline+str_iters+".csv"
    log_file = "Log/X01_Baseline_" + randOrline + str_iters + ".log"
    rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Baseline",html_file,cols=2)
    openLog(log_file,randOrline + str(iters))

    fake_geos =['RandA','RandB']
    #bins = [10,20,25,50,40,60,80,100,150,200,500,1000] #the last ones happen to be the ones we use for the plot
    densities = [1]  # the last ones happen to be the ones we use for the plot
    samples = [100]#,500,1000,2000,5000,10000]#,20000,50000]
    #samples = [200,20000]
    log(log_file,'Create samples')
    dic_sample_bin = {}
    for sample in samples:
        dic_sample_bin[sample] = {}
        dic_fake_1 = {'RandA': [], 'RandB': []}
        #dic_sample = {}
        for i in range(0,sample):
            if randOrline == 'spot':
                dic_fake_1['RandA'].append(50)
                dic_fake_1['RandB'].append(50)
            elif randOrline == 'line':
                dic_fake_1['RandA'].append(i)
                dic_fake_1['RandB'].append(i)
            elif randOrline == 'vert':
                dic_fake_1['RandA'].append(50)
                dic_fake_1['RandB'].append(i)
            elif randOrline == 'triangle':
                dic_fake_1['RandA'].append(i + random.randint(0,i))
                dic_fake_1['RandB'].append(i + random.randint(0,i))

        df_sample = pd.DataFrame.from_dict(dic_fake_1)

        for density in densities:
            log(log_file, 'samples=' + str(sample) + ' bins=' + str(density))
            dic_sample_bin[sample][density] =wcm.WilliamsDivergenceMaker(df_sample,fake_geos,density=density,log=1,norm=normed_corr,pval_iters=iters)
    # create a df with the coefficients
    dic_fake = {'stat': [],'p_value':[],'bins':[],'samples':[]}
    for sample in samples:
        for density in densities:
            log(log_file, 'Creating divergence dataframe=' + str(sample) + ' bin=' + str(density))
            cm = dic_sample_bin[sample][density]
            print(cm)
            div = cm.getCorrelation(['RandA', 'RandB'])
            stat, pvalue, A, D, B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
            dic_fake['stat'].append(stat)
            dic_fake['p_value'].append(pvalue)
            dic_fake['bins'].append(cm.bins)
            dic_fake['samples'].append(sample)
            # operation to normalise
    df_fake = pd.DataFrame.from_dict(dic_fake).round(4)
    print(df_fake)
    df_fake.to_csv(csv_file,index=False)
    #add a seaborn lineplot for all bins
    log(log_file, 'Make line plots')
    rep.addLineComment('Compare coefficients')
    fig, ax = plt.subplots()
    sns.lineplot(data=dic_fake,x='bins',y='stat',hue='samples',palette='Paired')
    plt.title('Change in correlations all bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    #add a seaborn lineplot for only up to 100
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake[df_fake['bins']<=100],x='bins',y='stat',hue='samples',palette='Paired')
    plt.title('Change in correlations <=100 bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    #add a seaborn lineplot for all bins on p-value
    rep.addLineComment('Compare coefficients')
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake,x='bins',y='p_value',hue='samples',palette='Paired')
    plt.title('Change in stat test all bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    #add a seaborn lineplot for only up to 100
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake[df_fake['bins']<=100],x='bins',y='p_value',hue='samples',palette='Paired')
    plt.title('Change in values <=100 bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)

    #add a seaborn lineplot where there is NO sin or cos involved
    rep.addLineComment('Data in the dataframes')
    for sample in samples:
        rep.addDataFrame(df_fake[df_fake['samples']==sample],title='Sample size=' + str(sample))

    log(log_file, 'Make distribution plots')
    rep.addLineComment('Plots of distributions')
    rep.changeColNumber(5)
    #plot_bins = [5,10,20,50]
    #plot_samples = [200,10000]
    for den in densities:
        for psample in samples:
            log(log_file, 'Plots bins=' + str(den) + ' size=' + str(psample))
            cm = dic_sample_bin[psample][den]
            cm_data = cm.data
            print(cm_data)
            df_rand = cm.randomiseData(cm_data,['RandA', 'RandB'])
            div = cm.getCorrelation(['RandA', 'RandB'])
            stat,pvalue,A,D,B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
            mean,sd,hist = div.p_mean,div.p_std,div.p_hist
            maxV = max(np.max(A),np.max(B))
            rep.addPlot2d(cm_data, 'scatter', title=str(round(stat,3)) + ' orig, bins=' + str(den), geo_x='RandA', geo_y='RandB', hue='RandA')
            rep.addPlot2d(df_rand, 'scatter',title='rand, size=' + str(psample), geo_x='RandA', geo_y='RandB', hue='RandA')
            #if len(hist['divergence'])>0:
            #    crit_val = round(cm.getCriticalValue('RandA','RandB',0.95),3)
            #    rep.addPlot1d(hist,'histogram',geo_x='divergence',title='mean=' + str(round(mean,3)) + ' sd=' + str(round(sd,3)) + ' crit5%=' + str(crit_val))
            #else:
            #    rep.addBoxComment(('No histogram calculated'))
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues',colourbar=False)
            rep.addSurface(D, 'Difference Data stat=' + str(round(stat,3)) + ' pvalue=' + str(round(pvalue,3)), cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds',colourbar=False)

    log(log_file, 'Finally print out')
    rep.printReport()


if __name__ == '__main__':
    globals()['randomBaseline'](sys.argv[1],sys.argv[2])

