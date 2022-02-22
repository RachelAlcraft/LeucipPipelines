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
    csv_correlations = "Csv/10_DivCorr_" + tag + ".csv"
    csv_summary = "Csv/10_SumCorr_" + tag + ".csv"
    log_file = "Log/10_OrderCorrellations_" + tag + ".log"
    openLog(log_file, csv_correlations)
    ###############################################################################################
    data = pd.read_csv(csv_final)
    geos = []
    for col in data.columns:
        #if ':' in col and 'info' not in col and '(' not in col and '{' not in col:
        #if ':' in col and 'info' not in col and '{' not in col:
        if ':' in col and 'info' not in col:
            geos.append(col)

    #geos = ['N:O','N:CA','CA:C','N:CA:C:N+1']

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        data = data.query('occupancy == 1')
        data = data.query('bfactor <= 10')
        geos_to_abs = ['CA:C:O:N+1', 'CA-1:C-1:N:CA', 'CA:C:N+1:CA+1']
        for gabs in geos_to_abs:
            data[gabs] = abs(data[gabs])
        aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
        data = data[data['aa'].isin(aa_list)]
        data = data[data['aa-1'].isin(aa_list)]
        data = data[data['aa+1'].isin(aa_list)]

    log(log_file,'Create Williams Coefficient Maker')
    wcc = wcm.WilliamsDivergenceMaker(data,geos,density=density,log=1,norm=True,pval_iters=0,delay_load=False)

    recreate_divergence = True
    if recreate_divergence:
        complete = wcc.getCoefficientsDataFrame()
        complete = complete.sort_values(by='stat', ascending=False)
        complete.to_csv(csv_correlations,index=False)
        print(complete)

        summary = complete.groupby(by='geoA').sum()
        summary = summary.sort_values(by='stat')
        summary.to_csv(csv_summary)

###################################################################################
orderCorrelations('High_GLY')
orderCorrelations('SYN_GLY')
orderCorrelations('High')
orderCorrelations('Redo')
orderCorrelations('Redo_GLY')




