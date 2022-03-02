'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
density = 5
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help

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
    csv_ideal = "Csv/PW_" + tag + "_IDEAL_01_Geometry.csv"
    csv_correlations = "Csv/10_DivCorr_" + tag + ".csv"
    csv_summary = "Csv/10_SumCorr_" + tag + ".csv"
    log_file = "Log/10_OrderCorrelations_" + tag + ".log"
    openLog(log_file, csv_correlations)
    ###############################################################################################
    data = pd.read_csv(csv_final)
    data_ideal = pd.read_csv(csv_ideal)
    geos = []
    for col in data.columns:
        if ':' in col and 'info' not in col and '(' not in col and '{' not in col:
        #if ':' in col and 'info' not in col and '{' not in col:
        #if ':' in col and 'info' not in col:
            geos.append(col)

    #geos = ['N:O','N:CA','CA:C','N:CA:C:N+1']

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        data = data.query('occupancy == 1')
        data = data.query('bfactor <= 10')
        data = data.query('`CA:C:N+1:CA+1` >= 100')
        data = data.query('`CA-1:C-1:N:CA` >= 100')
        geos_to_abs = ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1']
        for gabs in geos_to_abs:
            data[gabs] = abs(data[gabs])
        aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
        data = data[data['aa'].isin(aa_list)]
        data = data[data['aa-1'].isin(aa_list)]
        data = data[data['aa+1'].isin(aa_list)]

        for gabs in geos_to_abs:
            data_ideal[gabs] = abs(data_ideal[gabs])
        data_ideal = data_ideal.query('`CA:C:N+1:CA+1` >= 100')
        data_ideal = data_ideal.query('`CA-1:C-1:N:CA` >= 100')


    log(log_file,'Create Williams Coefficient Maker')
    wcc = wcm.WilliamsDivergenceMaker(data,geos,density=density,log=1,norm=True,pval_iters=0,delay_load=False)
    #wcc_ideal = wcm.WilliamsDivergenceMaker(data_ideal, geos, density=density, log=1, norm=True, pval_iters=0, delay_load=False)

    recreate_divergence = True
    if recreate_divergence:
        complete = wcc.getCoefficientsDataFrame()
        complete = complete.sort_values(by='stat', ascending=False)
        ideals = []
        for i in range(len(complete.index)):
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            tpl = wcc.compareTwoDistributions(data[[geoA,geoB]],data_ideal[[geoA,geoB]],geoA,geoB)
            ideals.append(tpl[0])
        complete['stat_ideal'] = ideals
        complete.to_csv(csv_correlations,index=False)
        print(complete)

        summary = complete.groupby(by='geoA').sum()
        summary = summary.sort_values(by='stat')
        summary.to_csv(csv_summary)

###################################################################################
#orderCorrelations('High_GLY')
#orderCorrelations('SYN_GLY')
#orderCorrelations('High')
#orderCorrelations('Redo')
orderCorrelations('High_GLY')
#orderCorrelations('Redo_GLY')
#orderCorrelations('Redo_GLY_IDEAL')
#orderCorrelations('High_GLY_IDEAL')




