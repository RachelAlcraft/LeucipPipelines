

import math
import random
import sys

import scipy
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import DataFrameMaker as dfm
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

def randomBaselineOne(str_iters,str_density,str_variance):
    normed_corr = False
    iters = int(str_iters)
    density = float(str_density)
    randOrline = 'three'#'line' #rand or line or covar
    vari = int(str_variance)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    html_file = "Html/01a_Baseline_"+randOrline+str_iters+'_'+str_density+'_'+str_variance+".html"
    csv_file = "Csv/01a_Baseline_"+randOrline+str_iters+'_'+str_density+'_'+str_variance+".csv"
    log_file = "Log/01a_Baseline_" + randOrline + str_iters+'_'+str_density+'_'+str_variance+".log"
    rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Density=" + str_density + ' variance=' + str_variance,html_file,cols=2)
    openLog(log_file,randOrline + str(iters))
    fake_geos =['RandA1','RandA2','RandB1','RandB2','RandC1','RandC2']
    samples = [200,500,750,1000,1500,2000,5000,10000]#,20000,50000]
    log(log_file,'Create samples')
    dic_sample_bin = {}
    for sample in samples:
        dic_sample_bin[sample] = {}
        #dic_fake_all = {'RandA1': [], 'RandA2': [],'RandB1': [], 'RandB2': [],'RandC1': [], 'RandC2': [],'RandD1': [],'RandD2': [], 'RandE1': [], 'RandE2': []}
        dic_fake_all = {'RandA1': [], 'RandA2': [], 'RandB1': [], 'RandB2': [], 'RandD1': [], 'RandD2': []}
        #dic_sample = {}
        for i in range(0,sample):
            #if randOrline == 'rand':
            dic_fake_all['RandA1'].append(random.randint(0,sample))
            dic_fake_all['RandA2'].append(random.randint(0,sample))
            #elif randOrline == 'line':
            dic_fake_all['RandB1'].append(i)
            dic_fake_all['RandB2'].append(i)
            #elif randOrline == 'triangle':
            #dic_fake_all['RandC1'].append(i + random.randint(0,i))
            #dic_fake_all['RandC2'].append(i + random.randint(0,i))
            #covar1
            dic_fake_all['RandD1'].append(i + random.randint(0, int(sample/vari)))
            dic_fake_all['RandD2'].append(i + random.randint(0, int(sample/vari)))
            #covar2
            #dic_fake_all['RandE1'].append(i + random.randint(0, sample/2))
            #dic_fake_all['RandE2'].append(i + random.randint(0, sample/2))

        df_sample = pd.DataFrame.from_dict(dic_fake_all)
        bin = int(math.sqrt(sample/density))
        log(log_file, 'samples=' + str(sample) + ' bins=' + str(bin))
        dic_sample_bin[sample] =wcm.WilliamsDivergenceMaker(df_sample,fake_geos,density=density,log=1,norm=normed_corr,pval_iters=iters,delay_load=True)
    # create a df with the coefficients
    dic_fake = {'stat': [],'p_value':[],'bins':[],'samples':[],'set':[],'density':[],'stat2':[]}
    for sample in samples:
        bin = int(math.sqrt(sample / density))
        log(log_file, 'Creating divergence dataframe=' + str(sample) + ' bin=' + str(bin))
        cm = dic_sample_bin[sample]
        divA = cm.getCorrelation(['RandA1', 'RandA2'])
        divB = cm.getCorrelation(['RandB1', 'RandB2'])
        #divC = cm.getCorrelation(['RandC1', 'RandC2'])
        divD = cm.getCorrelation(['RandD1', 'RandD2'])
        #divE = cm.getCorrelation(['RandE1', 'RandE2'])
        for div,nm in [[divA,'Random'],[divB,'Line'],[divD,'Var']]:
            stat, pvalue, A, D, B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
            density = sample/(bin**2)
            stat2 = math.log(stat*density)
            dic_fake['stat'].append(stat)
            dic_fake['stat2'].append(stat2)
            dic_fake['density'].append(density)
            dic_fake['p_value'].append(pvalue)
            dic_fake['bins'].append(bin)
            dic_fake['samples'].append(sample)
            dic_fake['set'].append(nm+str_density+str_variance)
            # operation to normalise
    df_fake = pd.DataFrame.from_dict(dic_fake).round(4)
    df_fake.to_csv(csv_file,index=False)
    #add a seaborn lineplot for all bins
    log(log_file, 'Make line plots')
    rep.addLineComment('Compare coefficients')
    rep.changeColNumber(3)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=dic_fake,x='samples',y='stat',hue='set',palette='Paired')
    plt.title('Change in metric as sample size varies')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=dic_fake, x='samples', y='stat2', hue='set', palette='Paired')
    plt.title('log(stat/density)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake,x='samples',y='p_value',hue='set',palette='Paired')
    plt.title('Change in stat test all bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake[df_fake['samples']<=2000], x='samples', y='stat', hue='set', palette='Paired')
    plt.title('Change in metric as sample size varies')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake[df_fake['samples']<=2000], x='samples', y='stat2', hue='set', palette='Paired')
    plt.title('log(stat/density)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=df_fake[df_fake['samples']<=2000], x='samples', y='p_value', hue='set', palette='Paired')
    plt.title('Change in stat test all bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################

    rep.changeColNumber(3)
    rep.addLineComment('Data in the dataframes')
    for sample in samples:
        rep.addDataFrame(df_fake[df_fake['samples']==sample],title='Sample size=' + str(sample))

    log(log_file, 'Make distribution plots')
    rep.addLineComment('Plots of distributions')
    rep.changeColNumber(6)

    for geoA,geoB in [['RandA1','RandA2'],['RandB1','RandB2'],['RandD1','RandD2']]:
        rep.addLineComment(geoA + geoB)
        for psample in samples:
            pbin = int(math.sqrt(sample / density))
            log(log_file, 'Plots bins=' + str(pbin) + ' size=' + str(psample))
            cm = dic_sample_bin[psample]
            cm_data = cm.data
            df_rand = cm.randomiseData(cm_data,[geoA, geoB])
            div = cm.getCorrelation([geoA, geoB])
            stat,pvalue,A,D,B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
            mean,sd,hist = div.p_mean,div.p_std,div.p_hist
            maxV = max(np.max(A),np.max(B))
            rep.addPlot2d(cm_data, 'scatter', title=str(round(stat,3)) + ' orig, bins=' + str(pbin), geo_x=geoA, geo_y=geoB, hue=geoA)
            rep.addPlot2d(df_rand, 'scatter',title='rand, size=' + str(psample), geo_x=geoA, geo_y=geoB, hue=geoA)
            if len(hist['divergence'])>0:
                crit_val = round(cm.getCriticalValue(geoA,geoB,0.95),3)
                rep.addPlot1d(hist,'histogram',geo_x='divergence',title='mean=' + str(round(mean,3)) + ' sd=' + str(round(sd,3)) + ' crit5%=' + str(crit_val))
            else:
                rep.addBoxComment(('No histogram calculated'))
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues',colourbar=False)
            rep.addSurface(D, 'Difference Data stat=' + str(round(stat,3)) + ' pvalue=' + str(round(pvalue,3)), cmin=-1*maxV, cmax=maxV, palette='RdBu',colourbar=False)
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds',colourbar=False)

    log(log_file, 'Finally print out')
    rep.printReport()


if __name__ == '__main__':
    globals()['randomBaselineOne'](sys.argv[1],sys.argv[2],sys.argv[3])

