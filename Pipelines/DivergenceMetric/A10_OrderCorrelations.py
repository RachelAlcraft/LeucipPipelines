'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help
import A0_Globals as glob

from LeucipPy import HtmlReportMaker as hrm
from LeucipPy import WilliamsDivergenceMaker as wcm

import pandas as pd
from datetime import datetime as dt

########################################################################################
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
########################################################################################
def orderCorrelations(tag):
    csv_final = "Csv/PW_" + tag + "_01_Geometry.csv"
    csv_ideal_paired = "Csv/PW_" + tag + "_IDEAL_PAIRED_01_Geometry.csv"
    csv_ideal_unpaired = "Csv/PW_" + tag + "_IDEAL_UNPAIRED_01_Geometry.csv"
    csv_correlations = "Csv/10_DivCorr_" + tag + ".csv"
    csv_summary = "Csv/10_SumCorr_" + tag + ".csv"
    log_file = "Log/10_OrderCorrelations_" + tag + ".log"
    openLog(log_file, csv_correlations)
    ###############################################################################################
    data = pd.read_csv(csv_final)
    data_ideal_paired = pd.read_csv(csv_ideal_paired)
    data_ideal_unpaired = pd.read_csv(csv_ideal_unpaired)
    geos = []
    for col in data.columns:
        if ':' in col and 'info' not in col and '(' not in col and '{' not in col and '4' not in col:
        #if ':' in col and 'info' not in col and '{' not in col:
        #if ':' in col and 'info' not in col:
            geos.append(col)

    #geos = ['N:O','N:CA','CA:C','N:CA:C:N+1']

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        data = glob.trimData(data,15,geos)
        data_ideal_paired = glob.trimData(data_ideal_paired, 15, geos)
        data_ideal_unpaired = glob.trimData(data_ideal_unpaired, 15, geos)
        geos_to_abs = glob.getGeosToAbs()
        for gabs in geos_to_abs:
            data[gabs] = abs(data[gabs])
            data_ideal_paired[gabs] = abs(data_ideal_paired[gabs])
            data_ideal_unpaired[gabs] = abs(data_ideal_unpaired[gabs])
        #data = data.query('`CA:C:N+1:CA+1` >= 100')
        #data = data.query('`CA-1:C-1:N:CA` >= 100')

    log(log_file,'Create Williams Coefficient Maker')
    wcc = wcm.WilliamsDivergenceMaker(data,geos,bins=10,log=1,norm=True,pval_iters=0,delay_load=False)
    print('Continuing...')
    #wcc_ideal = wcm.WilliamsDivergenceMaker(data_ideal, geos, density=density, log=1, norm=True, pval_iters=0, delay_load=False)

    recreate_divergence = True
    if recreate_divergence:
        complete = wcc.getCoefficientsDataFrame()
        complete = complete.sort_values(by='stat', ascending=False)
        ideals_p = []
        ideals_u = []
        ideals_up = []
        for i in range(len(complete.index)):
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            #print(geoA,geoB)
            tpla = wcc.compareTwoDistributions(data[[geoA,geoB]],data_ideal_paired[[geoA,geoB]],geoA,geoB)
            tplb = wcc.compareTwoDistributions(data[[geoA, geoB]], data_ideal_unpaired[[geoA, geoB]], geoA, geoB)
            tplc = wcc.compareTwoDistributions(data_ideal_paired[[geoA, geoB]], data_ideal_unpaired[[geoA, geoB]], geoA, geoB)
            ideals_p.append(tpla[0])
            ideals_u.append(tplb[0])
            ideals_up.append(tplc[0])
        complete['stat_ideal_paired'] = ideals_p
        complete['stat_ideal_unpaired'] = ideals_u
        complete['stat_paired_unpaired'] = ideals_up
        complete.to_csv(csv_correlations,index=False)
        print(complete)

        summary = complete.groupby(by='geoA').sum()
        summary = summary.sort_values(by='stat')
        summary.to_csv(csv_summary)

###################################################################################
#orderCorrelations('High_GLY')
#orderCorrelations('SYN_GLY')
orderCorrelations('High')
#orderCorrelations('Redo')
orderCorrelations('High_GLY')
#orderCorrelations('Redo_GLY')
#orderCorrelations('Redo_GLY_IDEAL')
#orderCorrelations('High_GLY_IDEAL')




