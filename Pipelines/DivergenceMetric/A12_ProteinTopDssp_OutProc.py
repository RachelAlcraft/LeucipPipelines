'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#csv_final = "C:/Dev/Github/LeucipPipelines/Pipelines/Geometry/04Compare/Csv/PW_High_GLY_02_Geometry.csv"
density = 5
num_top = 40
num_bottom = 10
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sys
sys.path.append('../1Library')
import Helpers as help

from LeucipPy import HtmlReportMaker as hrm

import pandas as pd
import numpy as np
import seaborn as sns
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
def proteinTopDssp(tag):
    csv_final = "Csv/PW_" + tag + "_01_Geometry.csv_dssp.csv"
    csv_correlations = "Csv/11_DivCorr_" + tag + ".csv"
    html_filename = 'Html/12_ProteinTop_dssp_' + tag + '.html'
    log_file = "Log/12_ProteinTop_dssp_" + tag + ".log"
    title = 'Williams Divergence from Trivial: Most and Least Correlated ' + tag
    openLog(log_file, csv_final)
    ###############################################################################################
    data = pd.read_csv(csv_final)
    complete = pd.read_csv(csv_correlations)
    complete = complete.sort_values(by='stat', ascending=False)
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
            data.loc[:,gabs] = abs(data.loc[:,gabs])
        aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
        data = data[data['aa'].isin(aa_list)]
        data = data[data['aa-1'].isin(aa_list)]
        data = data[data['aa+1'].isin(aa_list)]

    recreate_html = True
    if recreate_html:
        used_geos = []
        # create a custom pallette for dssp
        pal = sns.color_palette("bright")
        pal2 = []
        pal2.append('Gainsboro')
        for pl in pal:
            pal2.append(pl)
        customPalette = sns.set_palette(sns.color_palette(pal2))

        log(log_file, '### Creating html reports')
        rep_mak = hrm.HtmlReportMaker(title, html_filename, cols=4)
        rep_mak.addLineComment('The Ramachandran Plot')
        rep_mak.addPlot2d(data, 'seaborn', title='', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='aa-1', palette='tab20')
        rep_mak.addPlot2d(data, 'seaborn', title='', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='aa', palette='tab20')
        rep_mak.addPlot2d(data, 'seaborn', title='', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='aa+1', palette='tab20')
        rep_mak.addPlot2d(data, 'seaborn', title='', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='ss', palette=customPalette)

        rep_mak.addLineComment('Most correlated geos')
        i = -1
        while len(used_geos) < num_top:
            i += 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            stat = complete['stat'].values[i]
            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos:
                log(log_file,str(len(used_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB + ' stat=' + str(round(stat,3)))
                used_geos.append(geoA + '_' + geoB)

                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='aa-1',palette='tab20')
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='aa', palette='tab20')
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='aa+1', palette='tab20')
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='ss', palette=customPalette)


        log(log_file,'############# Least correlated ###########################')
        rep_mak.addLineComment('Least correlated geos')
        i = len(complete.index)
        last_geos = []
        while len(last_geos) < num_bottom:
            i -= 1
            geoA = complete['geoA'].values[i]
            geoB = complete['geoB'].values[i]
            stat = complete['stat'].values[i]

            if geoA + '_' + geoB not in used_geos and geoB + '_' + geoA not in used_geos:
                log(log_file,str(len(last_geos)) + ' ' + geoA + ' ' + geoB + '.........')
                rep_mak.addLineComment(geoA + ' | ' + geoB + ' stat=' + str(round(stat,3)))
                last_geos.append(geoA + '_' + geoB)

                #data = data.sort_values(by='aa-1', ascending=True)
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='aa-1', palette='tab20')
                #data = data.sort_values(by='aa', ascending=True)
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='aa', palette='tab20')
                #data = data.sort_values(by='aa+1', ascending=True)
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='aa+1', palette='tab20')
                #data = data.sort_values(by='ss', ascending=True)
                rep_mak.addPlot2d(data, 'seaborn', title='', geo_x=geoA, geo_y=geoB, hue='ss', palette=customPalette)

        log(log_file, 'Finally print out to ' + html_filename)
        rep_mak.printReport()

###################################################################################
if __name__ == '__main__':
    globals()['proteinTopDssp'](sys.argv[1])




