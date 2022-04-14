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
def orderCorrelations(set,tag,first,second,method,geobases,dimensions):
    csv_final = "Csv/PW_" + set + "_01_Geometry.csv"
    csv_ideal_paired = "Csv/PW_" + set + "_IDEAL_PAIRED_01_Geometry.csv"
    csv_ideal_unpaired = "Csv/PW_" + set + "_IDEAL_UNPAIRED_01_Geometry.csv"
    csv_correlations = "Csv/PW_" + set + '_' + first + '_'+second + '_'+tag + '_' + str(dimensions) + "_Metrics.csv"
    log_file = "Log/10_OrderCorrelations3d_" + set + tag + ".log"
    openLog(log_file, csv_correlations)
    ###############################################################################################
    data = pd.read_csv(csv_final)
    data_ideal_paired = pd.read_csv(csv_ideal_paired)
    data_ideal_unpaired = pd.read_csv(csv_ideal_unpaired)

    modify_csv = True
    if modify_csv:
        log(log_file, '### Applying filters to csv data')
        data = glob.trimDihs(data, 15)
        data_ideal_paired = glob.trimDihs(data_ideal_paired, 15)
        data_ideal_unpaired = glob.trimDihs(data_ideal_unpaired, 15)

    log(log_file, 'Create Williams Coefficient Maker')
    if first=='obs' and second=='none' and method=='diff':
        aw = awa.AlcraftWilliamsAssociation(data, bins=10, piters=0,loglevel=2)
    elif first=='obs' and second=='paired' and method=='diff':
        aw = awa.AlcraftWilliamsAssociation(data, data_ideal_paired, bins=10, piters=0,loglevel=2)
    elif first == 'obs' and second == 'unpaired' and method=='diff':
        aw = awa.AlcraftWilliamsAssociation(data, data_ideal_unpaired, bins=10, piters=0,loglevel=2)
    elif first == 'paired' and second == 'unpaired' and method=='diff':
        aw = awa.AlcraftWilliamsAssociation(data_ideal_paired, data_ideal_unpaired, bins=10, piters=0,loglevel=2)
    elif first == 'obs' and second == 'none' and method=='k-l':
        aw = awa.AlcraftWilliamsAssociation(data, bins=10, piters=0, method='k-l',loglevel=2)
    elif first == 'obs' and second == 'paired' and method=='k-l':
        aw = awa.AlcraftWilliamsAssociation(data, data_ideal_paired, bins=10, piters=0, method='k-l',loglevel=2)
    elif first == 'obs' and second == 'unpaired' and method=='k-l':
        aw = awa.AlcraftWilliamsAssociation(data, data_ideal_unpaired, bins=10, piters=0, method='k-l',loglevel=2)
    elif first == 'paired' and second == 'unpaired' and method=='k-l':
        aw = awa.AlcraftWilliamsAssociation(data_ideal_paired, data_ideal_unpaired, bins=10, piters=0, method='k-l',loglevel=2)

    recreate_divergence = True
    if recreate_divergence:
        dic_stats = {}
        geos = glob.getGeos(False, True)
        df_triv = aw.getStrongestAssociations(geobases,geos,dimensions,fraction=1/dimensions)
        df_triv.to_csv(csv_correlations,index=False)







###################################################################################
#orderCorrelations('High','obs','none','diff','N:CA:C',3)
#orderCorrelations('High','CO3','obs','none','diff',['C:O'],2)
#orderCorrelations('High','Search2','obs','none','diff',[],2)
#orderCorrelations('High','TAU3','obs','none','diff',['N:CA:C'],3)
#orderCorrelations('High','RAMA3','obs','none','diff',['C-1:N:CA:C','N:CA:C:N+1'],1)
orderCorrelations('High','Search','obs','none','diff',[],3)





