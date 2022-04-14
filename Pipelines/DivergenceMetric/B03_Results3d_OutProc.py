'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
num_top = 30
num_bottom = 5
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help
import A0_Globals as glob

from LeucipPy import HtmlReportMaker as hrm
from nDimAssociations import AlcraftWilliamsAssociation as awa
from nDimAssociations import ReportExport as re

import pandas as pd
import numpy as np
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
def proteinTop(tag,str_iters,stat):
    csv_final = "Csv/PW_" + tag + "_01_Geometry.csv"
    csv_ideal_paired = "Csv/PW_" + tag + "_IDEAL_PAIRED_01_Geometry.csv"
    csv_ideal_unpaired = "Csv/PW_" + tag + "_IDEAL_UNPAIRED_01_Geometry.csv"
    csv_correlations = "Csv/10_DivCorr3d_" + tag + ".csv"
    print(csv_correlations)
    html_filename = 'Html/B03_ProteinTop3d_' + tag + str_iters + stat + '.html'
    log_file = "Log/B03_ProteinTop3d_" + tag + str_iters + stat+ ".log"
    title = 'Alcraft-Williams Associations: 3d Associations ' + tag
    openLog(log_file, str_iters)
    ###############################################################################################
    csvA = csv_final
    csvB = ''
    method = 'diff'
    if stat == 'stat_obs_paired':
        csvB = csv_ideal_paired
    elif stat == 'stat_obs_unpaired':
        csvB = csv_ideal_unpaired
    elif stat == 'stat_paired_unpaired':
        csvA = csv_ideal_paired
        csvB = csv_ideal_unpaired
    elif stat == 'kl_trivial':
        method = 'k-l'
    elif stat == 'kl_obs_paired':
        csvB = csv_ideal_paired
        method = 'k-l'
    elif stat == 'kl_obs_unpaired':
        csvB = csv_ideal_unpaired
        method = 'k-l'
    elif stat == 'kl_paired_unpaired':
        csvA = csv_ideal_paired
        csvB = csv_ideal_unpaired
        method = 'k-l'

    iters = int(str_iters)
    dataA = pd.read_csv(csvA)
    dataB = None
    if csvB != '':
        dataB = pd.read_csv(csvB)

    modify_csv = True
    if modify_csv:
        log(log_file,'### Applying filters to csv data')
        dataA = glob.trimDihs(dataA, 15)
        if csvB != '':
            dataB = glob.trimDihs(dataB, 15)

    log(log_file,'Create Williams Coefficient Maker')
    aw = None
    if csvB != '':
        aw = awa.AlcraftWilliamsAssociation(dataA,dataB, piters=iters,method=method)
    else:
        aw = awa.AlcraftWilliamsAssociation(dataA,piters=iters,method=method)
    complete = pd.read_csv(csv_correlations)
    complete = complete.sort_values(by=stat, ascending=False)
    recreate_html = True
    dont_use = ['C-1:CB-1','CA-1:C-1','N-1:C-1','N-1:CA-1:C-1','CB:CA:C','N-1:CA-1','N-1:CB-1','CA-1:CB-1']

    if recreate_html:
        used_geos = []
        log(log_file,'### Creating html reports')
        rep_mak = re.ReportExport(title, html_filename, cols=4)
        rep_mak.addLineComment('Most correlated geos')
        i = -1
        count_up = 0
        while count_up < num_top and i < len(complete.index):
            i += 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            geoC = complete['geoC'].values[i]
            # all id tags
            tag1 = geoA + '_' + geoB + '_' + geoC
            tag2 = geoA + '_' + geoC + '_' + geoB
            tag3 = geoB + '_' + geoA + '_' + geoC
            tag4 = geoB + '_' + geoC + '_' + geoA
            tag5 = geoC + '_' + geoA + '_' + geoB
            tag6 = geoC + '_' + geoB + '_' + geoA

            if tag1 not in used_geos and geoA not in dont_use and geoB not in dont_use and geoC not in dont_use:
                count_up += 1
                log(log_file,str(count_up) + ' ' + geoA + ' ' + geoB + ' ' + geoC + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB + '|' + geoC)
                used_geos.append(tag1)
                used_geos.append(tag2)
                used_geos.append(tag3)
                used_geos.append(tag4)
                used_geos.append(tag5)
                used_geos.append(tag6)
                div = aw.addAssociation([geoA,geoB,geoC])
                stat, pvalue, A, D, B = div.metric, div.pvalue, div.matA, div.matDiff, div.matB
                info = geoA + ' | ' + geoB + ' | ' + geoC + ' , Metric = ' + str(round(stat,4)) + " pvalue = " + str(round(pvalue,4))
                rep_mak.addLineComment(info)
                rep_mak.addPlot2d(dataA, 'scatter', title='Compare', geo_x=geoA, geo_y=geoB, hue=geoC)
                rep_mak.addPlot2d(dataA, 'scatter', title='Compare', geo_x=geoB, geo_y=geoC, hue=geoA)
                rep_mak.addPlot2d(dataA, 'scatter', title='Compare', geo_x=geoC, geo_y=geoA, hue=geoB)


        log(log_file,'############# Least correlated ###########################')
        rep_mak.addLineComment('Least correlated geos')
        i = len(complete.index)
        count_down = 0
        while count_down < num_bottom and i >= 0:
            i -= 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            geoC = complete['geoC'].values[i]

            tag1 = geoA + '_' + geoB + '_' + geoC
            tag2 = geoA + '_' + geoC + '_' + geoB
            tag3 = geoB + '_' + geoA + '_' + geoC
            tag4 = geoB + '_' + geoC + '_' + geoA
            tag5 = geoC + '_' + geoA + '_' + geoB
            tag6 = geoC + '_' + geoB + '_' + geoA

            if tag1 not in used_geos and geoA not in dont_use and geoB not in dont_use and geoC not in dont_use:
                count_down += 1
                log(log_file,str(len(used_geos)) + ' ' + geoA + ' ' + geoB + ' ' + geoC + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB + '|' + geoC)
                used_geos.append(tag1)
                used_geos.append(tag2)
                used_geos.append(tag3)
                used_geos.append(tag4)
                used_geos.append(tag5)
                used_geos.append(tag6)

                log(log_file, str(count_down) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB)
                div = aw.addAssociation([geoA, geoB,geoC])
                stat, pvalue, A, D, B = div.metric, div.pvalue, div.matA, div.matDiff, div.matB
                info = geoA + ' | ' + geoB + ' | ' + geoC + ' , Metric = ' + str(round(stat,4)) + " pvalue = " + str(round(pvalue,4))
                rep_mak.addLineComment(info)
                print(dataA)
                rep_mak.addPlot2d(dataA, 'scatter', title='Compare', geo_x=geoA, geo_y=geoB, hue=geoC)
                rep_mak.addPlot2d(dataA, 'scatter', title='Compare', geo_x=geoB, geo_y=geoC, hue=geoA)
                rep_mak.addPlot2d(dataA, 'scatter', title='Compare', geo_x=geoC, geo_y=geoA, hue=geoB)

        log(log_file, 'Finally print out to ' + html_filename)
        rep_mak.printReport()

###################################################################################
if __name__ == '__main__':
    globals()['proteinTop'](sys.argv[1],sys.argv[2],sys.argv[3])




