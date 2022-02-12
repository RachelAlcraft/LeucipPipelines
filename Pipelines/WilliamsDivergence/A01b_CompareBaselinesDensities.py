

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

def run():
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    html_file = "Html/01b_CompareBaselines.html"
    log_file = "Log/01b_CompareBaselines.log"
    openLog(log_file, 'Compare baselines')
    csv_files = []
    csv_files.append("Csv/01a_Baseline_three500_0.5.csv")
    csv_files.append("Csv/01a_Baseline_three500_1.csv")
    csv_files.append("Csv/01a_Baseline_three500_5.csv")
    csv_files.append("Csv/01a_Baseline_three500_10.csv")
    csv_files.append("Csv/01a_Baseline_three500_20.csv")
    frames = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        frames.append(df)
    fake_df = pd.concat(frames)
    fake_df.index = range(0,len(fake_df.index))
    print(fake_df)

    rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Compare Baselines",html_file,cols=2)

    #add a seaborn lineplot for all bins
    log(log_file, 'Make line plots')
    rep.addLineComment('Compare coefficients')
    rep.changeColNumber(3)
    ############################################
    #rep.addPlot2d(fake_df,'seaborn',geo_x='samples',geo_y='stat',hue='set')
    hue_order = fake_df['set'].unique()
    hue_order.sort()
    colours = ["#F9EBEA", "#E6B0AA","#922B21", "#641E16","#CD6155","#F7DC6F", "#ABEBC6", "#1F618D", "#1B4F72", "#5499C7","#D4E6F1", "#ABEBC6", "#1E8449", "#145A32", "#52BE80"]

    fig, ax = plt.subplots()
    sns.lineplot(data=fake_df,x='samples',y='stat',hue='set',hue_order=hue_order,palette=colours)
    plt.title('Change in metric as sample size varies')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=fake_df, x='samples', y='stat2', hue='set',hue_order=hue_order,palette=colours)
    plt.title('log(stat/density)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=fake_df,x='samples',y='p_value',hue='set',hue_order=hue_order,palette=colours)
    plt.title('Change in stat test all bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=fake_df[fake_df['samples']<=2000], x='samples', y='stat', hue='set',hue_order=hue_order,palette=colours)
    plt.title('Change in metric as sample size varies')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=fake_df[fake_df['samples']<=2000], x='samples', y='stat2', hue='set',hue_order=hue_order,palette=colours)
    plt.title('log(stat/density)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################
    fig, ax = plt.subplots()
    sns.lineplot(data=fake_df[fake_df['samples']<=2000], x='samples', y='p_value', hue='set',hue_order=hue_order,palette=colours)
    plt.title('Change in stat test all bins')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.ylim(0, 1)
    rep.addPlotOnly(fig, ax)
    ############################################

    rep.changeColNumber(3)
    rep.addLineComment('Data in the dataframes')
    for sample in fake_df['samples'].unique():
        rep.addDataFrame(fake_df[fake_df['samples']==sample],title='Sample size=' + str(sample))


    log(log_file, 'Finally print out')
    rep.printReport()

############################################################################
run()