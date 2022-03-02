'''
Rachel Alcraft: 24/02/2022
Nice looking Ramachandran Plot
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
from Pipelines.Demonstrations import A0Functions

recreate_csv,modify_csv,recreate_html = False,True,True
# directory and file names
#directory of the pdb files
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data_redo/'
#full path, no path saves to current directory
csv_final = "csv/PrettyRamachandran.csv"
html_filename = 'Html/PrettyRamachandran.html'
#the geometric measures for geometry calculations
geos = ['N:CA:C:N+1','C-1:N:CA:C','C-1:C','N:N+1','N:O','N:CA:C']
title = 'Pretty Ramachandran Plot' # used at the header of the html report
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
###############################################################################################
if recreate_csv:
    print('### Load structures from BioPython #############')
    strucs = bpm.loadPdbStructures([],dir,extension='ent',prefix='pdb',log=2)
    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(geos, log=1)
    print('### Save dataframe ###############################')
    data.to_csv(csv_final , index=False)
    print("Saved to", csv_final)

data = pd.read_csv(csv_final)
if modify_csv:
    print('### Applying filters to csv data')

if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=2)
    print('--- data scatter')
    rep_mak.addPlot2d(data,'scatter',geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='N:CA:C', alpha=0.6, palette='jet_r')

    #https: // seaborn.pydata.org / examples / layered_bivariate_plot.html
    # Draw a combo histogram and scatterplot with density contours
    for df in [data]:
        print('--- scatters')
        x = df['C-1:N:CA:C']
        y = df['N:CA:C:N+1']



        f, ax = plt.subplots(figsize=(6, 6))
        sns.scatterplot(x=x, y=y, s=2, color=".15")
        sns.histplot(x=x, y=y, bins=200, pthresh=.1, cmap="Spectral_r")
        #print('--- skde 1')
        #xkde, ykde = A0Functions.make_Kde_able(x, y, 50, 10)
        sns.kdeplot(x=x, y=y, levels=10, color="AliceBlue", linewidths=0.5, bw_method=0.2)
        rep_mak.addPlotOnly(f,ax)

    rep_mak.printReport()


