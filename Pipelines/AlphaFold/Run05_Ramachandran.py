'''
Rachel Alcraft: 31/01/2022
Script to determine simply the bond and angle lengths
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
import numpy as np
from scipy.stats import gaussian_kde

from Pipelines.AlphaFold import A0Functions

recreate_csv,recreate_html = False,True
tag = 'Human'
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_alpha/' + tag + '/'
csv_output = "Csv/Ramachandran_" + tag + ".csv"
html_output = 'Html/' + tag + "_Ramachandran.html"
geos = ['N:CA','CA:C','C:O','C:N+1','N:CA:C','C-1:N:CA:C','N:CA:C:N+1']
title = 'Geometry Inspection, Alpha Fold' # used at the header of the html report
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
###############################################################################################
starts = []
count = 1000
start = 0
while start < 23391:
    starts.append(start)
    start += count

#starts = starts[-1:]
print(starts)

if recreate_csv:
    print('### Load structures from BioPython #############')
    # there are 23,391 humans
    for start in starts:
        strucs = bpm.loadPdbStructures([],dir,extension='pdb', prefix='',log=2, start=start, count=count)
        print('### Creating dataframe for correlations #############')
        geo_mak = dfm.GeometryMaker(strucs, log=1)
        data = geo_mak.calculateGeometry(geos, log=1)
        print('### Save dataframe ###############################')
        csv_output = "Csv/Ramachandran_" + tag + '_' + str(start) + ".csv"
        data.to_csv(csv_output , index=False)
        print("Saved to", csv_output)

else:
    print('### Load from csv file #############')
    csv_input = "Csv/Ramachandran_" + tag + '_' + str(0) + ".csv"
    data = pd.read_csv(csv_input)
    for start in starts[1:]:
        csv_input = "Csv/Ramachandran_" + tag + '_' + str(start) + ".csv"
        datanew = pd.read_csv(csv_input)
        data = data.append(datanew,ignore_index=True)
        print(len(data.index))
modify_csv = True
if modify_csv:
    print('### Applying filters to csv data')
    data['Probability'] = data['bfactor']
    dataHigh = data.query('Probability >= 90')
    dataMiddle = data.query('Probability >= 50')
    dataMiddle = dataMiddle.query('Probability < 90')
    dataLow = data.query('Probability < 50')


if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_output, cols=3)
    rep_mak.addBoxComment('High probability >=90%')
    rep_mak.addBoxComment('Medium probability 50%-90%')
    rep_mak.addBoxComment('Low probability <90%')
    print('--- dataHigh scatter')
    rep_mak.addPlot2d(dataHigh,'scatter',geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])
    print('--- dataMiddle scatter')
    rep_mak.addPlot2d(dataMiddle, 'scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])
    print('--- dataLow scatter')
    rep_mak.addPlot2d(dataLow, 'scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])

    #https: // seaborn.pydata.org / examples / layered_bivariate_plot.html
    # Draw a combo histogram and scatterplot with density contours
    for df in [dataHigh,dataMiddle,dataLow]:
        print('--- probability')
        x = df['C-1:N:CA:C']
        y = df['N:CA:C:N+1']


        f, ax = plt.subplots(figsize=(6, 6))
        sns.scatterplot(x=x, y=y, s=2, color=".15")
        sns.histplot(x=x, y=y, bins=200, pthresh=.1, cmap="Spectral_r")
        #print('--- skde 1')
        #xn, yn = A0Functions.make_Kde_able(x, y, 50, 10000)
        df_resampled = df.sample(frac=0.0001,replace=True)
        print('Original length=', len(df.index))
        print('Resampled length=',len(df_resampled.index))
        xn = df_resampled['C-1:N:CA:C']
        yn = df_resampled['N:CA:C:N+1']

        sns.kdeplot(x=xn, y=yn, levels=10, color="AliceBlue", linewidths=0.5,bw_method=0.2)
        rep_mak.addPlotOnly(f,ax)

    rep_mak.printReport()


