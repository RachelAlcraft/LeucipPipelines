'''
Rachel Alcraft: 31/01/2022
Script using LeucipPy Protein Geometry Library
This script:
 1) Loads BioPython structures from pdb files in a directory
 2) Creates a dataframe that correlates chosen geoemtric measures for all those structures, you can use this to process the data elsewhere
 3) Displays chosen correlations of those geos in an html report driven by matplotlib and seaborn

 The leucipPy library can be installed from TestPyPi with the command

 pip install -i https://test.pypi.org/simple/ LeucipPy-pkg-RachelAlcraft

 There is a google colab which demonstrates it here

 https://colab.research.google.com/drive/1Ut5HXQQgE3sNuAMo7IVR3aqwUVOYDy12#scrollTo=2ubfgnaNP_Gs
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
recreate_csv,modify_csv,recreate_html = False,False,True
# directory and file names
#directory of the pdb files
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data_redo/'
#full path, no path saves to current directory
csv_final = "Csv/GeometyInspection.csv"
html_filename = 'Html/GeometryInspection.html'
#the geometric measures for geometry calculations
geos = ['N:CA','CA:C','C:O','C:N+1','N:(N,S,O)+2','O:(N,S,O)+2']
title = 'Geometry Inspection' # used at the header of the html report
# Help on geos:
### distances, angles or dihedrals/improper angles by e.g. N:CA N:CA:C or N:CA:C:N+1
# Help on contact search:
### via brackets where {CB} means the nearest carbon beta atom and (C) nearest carbon atom of any type
### A list means the nearest og any in the list, e.g. O:(N,O) means the nearest of any N or O atoms
### You can specify not within x residues, eg O:(N,O)+2 means at least 2 residues away
### You can specify second closest eg SG:{SG@2} means second closest SG atom and O:{CB@2}+3 second clostes carbon beta not within 3 residues
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import DataFrameMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
import pandas as pd
###############################################################################################
if recreate_csv:
    print('### Load structures from BioPython #############')
    strucs = bpm.loadPdbStructures([],dir,extension='ent',prefix='pdb',log=2)
    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.DataFrameMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(geos, log=1)
    print('### Save dataframe ###############################')
    data.to_csv(csv_final , index=False)
    print("Saved to", csv_final)

data = pd.read_csv(csv_final)
if modify_csv:
    print('### Applying filters to csv data')
    # for alphafold we can rename the bfactor column  probability
    data['probability'] = data['bfactor']
    # it looks better when higher probability is sorted to the top
    data = data.sort_values(by='Probability', ascending=True)

if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=3)
    # Some scatter plots are added to make better sense of the data, the colour scheme is to approcimate the AlphaFold probabilty
    rep_mak.addLineComment('Scatter Plots')

    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'C:O-1', hue='O:CA-1', title='', palette='jet_r',xrange=[-180,180],yrange=[2.7,5.5])
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'C:CA-1', hue='O:CA-1', title='', palette='jet_r',xrange=[-180,180],yrange=[2.7,5.5])
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'O-1:N+1', hue='O:CA-1', title='', palette='jet_r',xrange=[-180,180],yrange=[2.7,5.5])

    rep_mak.addPlot2d(data,'scatter','C-1:N:CA:C','N:CA:C:N+1',hue='O:CA-1',title='',palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='O:O-1', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='O:CA-1', title='', palette='jet_r')

    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='C:O-1', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='C:CA-1', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='C:CA+1', title='', palette='jet_r')

    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='O-1:N+1', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='O-1:CA+1', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='O-1:C+1', title='', palette='jet_r')

    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='C-1:N:CA', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N:CA:C', title='', palette='jet_r')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='CA:C:N+1', title='', palette='jet_r')


    rep_mak.printReport()


