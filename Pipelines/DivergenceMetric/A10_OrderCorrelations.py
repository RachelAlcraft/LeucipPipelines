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

from nDimAssociations import AlcraftWilliamsAssociation as awa

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
    geos = glob.getGeos(False,True)

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        data = glob.trimDihs(data,15)
        data_ideal_paired = glob.trimDihs(data_ideal_paired, 15)
        data_ideal_unpaired = glob.trimDihs(data_ideal_unpaired,15)

    log(log_file,'Create Williams Coefficient Maker')
    aw_triv = awa.AlcraftWilliamsAssociation(data,bins=10,piters=0)
    aw_paired = awa.AlcraftWilliamsAssociation(data, data_ideal_paired,bins=10, piters=0)
    aw_unpaired = awa.AlcraftWilliamsAssociation(data, data_ideal_unpaired,bins=10, piters=0)
    aw_pairunpaired = awa.AlcraftWilliamsAssociation(data_ideal_paired, data_ideal_unpaired, bins=10, piters=0)
    kl_triv = awa.AlcraftWilliamsAssociation(data, bins=10, piters=0,method='k-l')
    kl_paired = awa.AlcraftWilliamsAssociation(data, data_ideal_paired, bins=10, piters=0,method='k-l')
    kl_unpaired = awa.AlcraftWilliamsAssociation(data, data_ideal_unpaired, bins=10, piters=0,method='k-l')
    kl_pairunpaired = awa.AlcraftWilliamsAssociation(data_ideal_paired, data_ideal_unpaired, bins=10, piters=0,method='k-l')
    print('Continuing...')

    recreate_divergence = True
    if recreate_divergence:
        dic_stats = {}
        dic_stats['geoA'] = []
        dic_stats['geoB'] = []
        dic_stats['stat_trivial'] = []
        dic_stats['stat_obs_paired'] = []
        dic_stats['stat_obs_unpaired'] = []
        dic_stats['stat_paired_unpaired'] = []
        dic_stats['kl_trivial'] = []
        dic_stats['kl_obs_paired'] = []
        dic_stats['kl_obs_unpaired'] = []
        dic_stats['kl_paired_unpaired'] = []
        for a in range(len(geos)):
            geoA = geos[a]
            for b in range(a+1,len(geos)):
                geoB = geos[b]
                if (a+b)%100 == 0:
                    print(a,b,'/',len(geos),geoA,geoB)
                s1 = aw_triv.addAssociation([geoA,geoB])
                s2 = aw_paired.addAssociation([geoA,geoB])
                s3 = aw_unpaired.addAssociation([geoA,geoB])
                s4 = aw_pairunpaired.addAssociation([geoA,geoB])
                s5 = kl_triv.addAssociation([geoA, geoB])
                s6 = kl_paired.addAssociation([geoA, geoB])
                s7 = kl_unpaired.addAssociation([geoA, geoB])
                s8 = kl_pairunpaired.addAssociation([geoA, geoB])
                dic_stats['stat_trivial'].append(s1.metric)
                dic_stats['stat_obs_paired'].append(s2.metric)
                dic_stats['stat_obs_unpaired'].append(s3.metric)
                dic_stats['stat_paired_unpaired'].append(s4.metric)
                dic_stats['kl_trivial'].append(s5.metric)
                dic_stats['kl_obs_paired'].append(s6.metric)
                dic_stats['kl_obs_unpaired'].append(s7.metric)
                dic_stats['kl_paired_unpaired'].append(s8.metric)
                dic_stats['geoA'].append(geoA)
                dic_stats['geoB'].append(geoB)

        complete = pd.DataFrame.from_dict(dic_stats)
        complete.to_csv(csv_correlations,index=False)
        summary = complete.groupby(by='geoA').sum()
        summary = summary.sort_values(by='stat_trivial')
        summary.to_csv(csv_summary)

###################################################################################
#orderCorrelations('High_GLY')
#orderCorrelations('SYN_GLY')
orderCorrelations('High')
#orderCorrelations('Redo')
#orderCorrelations('High_GLY')
#orderCorrelations('Redo_GLY')
#orderCorrelations('Redo_GLY_IDEAL')
#orderCorrelations('High_GLY_IDEAL')




