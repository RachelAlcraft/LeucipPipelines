

import math
import sys

from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsDivergenceMaker as wcm
import A0_Globals as globals
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

def run(str_density,str_bins,str_iters):
    normed_corr = True
    density = float(str_density)
    iters = int(str_iters)
    bins = int(str_bins)
    randOrline = 'three'#'line' #rand or line or covar
    varis = [0,0.1,0.5,1,10]#,50,500]
    resample = True
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    html_file = "Html/A01_Baseline_"+randOrline+str(iters)+'_'+str(density)+".html"
    csv_file = "Csv/A01_Baseline_"+randOrline+str(iters)+'_'+str(density)+".csv"
    log_file = "Log/A01_Baseline_" + randOrline + str(iters)+'_'+str(density)+".log"
    rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Variance Comparison, Density=" + str(density),html_file,cols=2)
    openLog(log_file,randOrline + str(iters))
    samples = [200,500,750,1000,1500,2000,5000,10000]
    fake_geos = []
    geo_pairs = []
    log(log_file,'Create samples')
    div_per_sample = {}

    for sample in samples:
        dic_fake_all = {}
        for vi in varis:
            log(log_file, 'samples=' + str(sample) + ' vi=' + str(vi))
            tag = str(vi) + '_'+ str(sample)
            fake_geos.append(tag+'A')
            fake_geos.append(tag+'B')
            geo_pairs.append([sample,vi,tag+'A',tag+'B'])
            dic_fake_all[tag+'A'] = []
            dic_fake_all[tag+'B'] = []
            for i in range(0,sample):
                l = i % 20
                count = 10  # samplesize
                dic_fake_all[tag+'A'].append(np.random.normal(l, int(count * vi)))
                dic_fake_all[tag+'B'].append(np.random.normal(l, int(count * vi)))
                #dic_fake_all[tag+'A'].append(np.random.normal(i, int(sample*vi)))
                #dic_fake_all[tag + 'B'].append(np.random.normal(i, int(sample * vi)))

        #print(dic_fake_all)
        df_sample = pd.DataFrame.from_dict(dic_fake_all)
        print(df_sample.columns)
        if str_density == '0':
            density = sample/(bins*bins)
        div_per_sample[sample] =wcm.WilliamsDivergenceMaker(df_sample,fake_geos,density=density,log=1,norm=normed_corr,pval_iters=iters,delay_load=True,p_resample=resample)
        print('######',sample)
        print(div_per_sample,sample)
        print('###############')
    # create a df with the coefficients
    dic_fake = {'stat': [],'p_value':[],'bins':[],'samples':[],'set':[],'density':[],'stat2':[],'random':[]}
    for sample,vi,geoA,geoB in geo_pairs:
        log(log_file, 'Creating divergence dataframe=' + str(sample) + geoA + ' ' + geoB)
        dm = div_per_sample[sample]
        div = dm.getCorrelation([geoA,geoB])
        stat, pvalue, A, D, B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
        stat2 = math.log(stat*density)
        dic_fake['stat'].append(stat)
        dic_fake['stat2'].append(stat2)
        dic_fake['density'].append(density)
        dic_fake['p_value'].append(pvalue)
        dic_fake['bins'].append(dm.bins)
        dic_fake['samples'].append(sample)
        dic_fake['random'].append(vi)
        dic_fake['set'].append(geoA)
                # operation to normalise
    df_fake = pd.DataFrame.from_dict(dic_fake).round(4)
    df_fake.to_csv(csv_file,index=False)
    #add a seaborn lineplot for all bins
    log(log_file, 'Make line plots')
    rep.addLineComment('Compare coefficients')
    rep.changeColNumber(2)

    #rep.addPlot2d(df_fake,'scatter',geo_x='samples',geo_y='stat',hue='random',yrange=[0,1])
    #rep.addPlot2d(df_fake, 'scatter', geo_x='random', geo_y='stat', hue='samples',yrange=[0,1])
    #rep.addPlot2d(df_fake[df_fake['samples']==1000], 'seaborn', geo_x='random', geo_y='stat', hue='set',yrange=[0,1])

    ############################################
    ys = ['stat','p_value']
    for y in ys:
        fig, ax = plt.subplots()
        sns.lineplot(data=df_fake,x='random',y=y,hue='samples',palette='tab10')
        plt.title('Change at different randomness, per sample size')
        plt.legend(title='sample size',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.ylim(0,1)
        rep.addPlotOnly(fig, ax)
        ############################################
        fig, ax = plt.subplots()
        sns.lineplot(data=df_fake, x='samples', y=y, hue='random',palette='tab10')
        plt.title('Change over sample size, per randomness')
        plt.legend(title='randomness',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.ylim(0, 1)
        rep.addPlotOnly(fig, ax)
        ##############################

    rep.changeColNumber(3)
    rep.addLineComment('Data in the dataframes')
    for sample in samples:
        rep.addDataFrame(df_fake[df_fake['samples']==sample],title='Sample size=' + str(sample))

    log(log_file, 'Make distribution plots')
    rep.addLineComment('Plots of distributions')

    rep.changeColNumber(7)
    for sample,vi, geoA, geoB in geo_pairs:
        #rep.addLineComment(geoA + geoB, + 'Samples = ')
        if sample == 1000:
            log(log_file, 'Plots ' + geoA + ' ' + geoB)
            dm = div_per_sample[sample]
            cm_data = dm.data
            df_rand = dm.randomiseData(cm_data,[geoA, geoB])
            df_samp = dm.resampleData(cm_data,[geoA,geoB])
            print(df_samp)
            div = dm.getCorrelation([geoA, geoB])
            stat,pvalue,A,D,B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
            mean,sd,hist = div.p_mean,div.p_std,div.p_hist
            maxV = max(np.max(A),np.max(B))
            rep.addPlot2d(cm_data, 'scatter', title=str(round(stat,3)) + ' orig vari=' + str(vi), geo_x=geoA, geo_y=geoB, hue=geoA)
            rep.addPlot2d(df_rand, 'scatter',title='rand, size=' + str(sample), geo_x=geoA, geo_y=geoB, hue=geoA)
            rep.addPlot2d(df_samp, 'scatter', title='resampled, size=' + str(sample), geo_x=geoA, geo_y=geoB, hue=geoA)
            if len(hist['divergence'])>0:
                crit_val = round(dm.getCriticalValue(geoA,geoB,0.95),3)
                rep.addPlot1d(hist,'histogram',geo_x='divergence',title='mean=' + str(round(mean,3)) + ' sd=' + str(round(sd,3)) + ' crit5%=' + str(crit_val))
            else:
                rep.addBoxComment(('No histogram calculated'))
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues',colourbar=False)
            rep.addSurface(D, 'Difference Data stat=' + str(round(stat,3)) + ' pvalue=' + str(round(pvalue,3)), cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds',colourbar=False)

    log(log_file, 'Finally print out')
    rep.printReport()


if __name__ == '__main__':
    globals()['run'](sys.argv[1],sys.argv[2],sys.argv[3])
