'''
Rachel Alcraft: 31/01/2022
Script to determine simply the bond and angle lengths
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
recreate_csv,modify_csv,recreate_html = True,True,True
dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data_redo/'
csv_output = "Csv/BondAngleLengths.csv"
html_output = 'Html/BondAngleLengths.html'
geos = ['N:CA','CA:C','C:O','C:N+1','N:CA:C','C-1:N:CA:C','N:CA:C:N+1']
title = 'Geometry Inspection' # used at the header of the html report
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
import pandas as pd
###############################################################################################
if recreate_csv:
    print('### Load structures from BioPython #############')
    strucs = bpm.loadPdbStructures([],dir,extension='ent',prefix='pdb',log=2)
    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(geos, log=1)
    print('### Save dataframe ###############################')
    data.to_csv(csv_output , index=False)
    print("Saved to", csv_output)

data = pd.read_csv(csv_output)
if modify_csv:
    print('### Applying filters to csv data')
    data = data.query('bfactor < 10')

if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_output, cols=4)
    rep_mak.addPlot1d(data, 'histogram', geo_x='N:CA',hue='pdb_code')
    rep_mak.addPlot1d(data, 'histogram', geo_x='CA:C', hue='pdb_code')
    rep_mak.addPlot1d(data, 'histogram', geo_x='C:O', hue='pdb_code')
    rep_mak.addPlot1d(data, 'histogram', geo_x='C:N+1', hue='pdb_code')
    rep_mak.addSeries(data['N:CA'].describe(), '', True)
    rep_mak.addSeries(data['CA:C'].describe(), '', True)
    rep_mak.addSeries(data['C:O'].describe(), '', True)
    rep_mak.addSeries(data['C:N+1'].describe(), '', True)

    rep_mak.addPlot1d(data, 'histogram', geo_x='N:CA:C', hue='pdb_code')
    rep_mak.addPlot1d(data, 'histogram', geo_x='C-1:N:CA:C', hue='pdb_code')
    rep_mak.addPlot1d(data, 'histogram', geo_x='N:CA:C:N+1', hue='pdb_code')
    rep_mak.addBoxComment('')
    rep_mak.addSeries(data['N:CA:C'].describe(), '', True)
    rep_mak.addSeries(data['C-1:N:CA:C'].describe(), '', True)
    rep_mak.addSeries(data['N:CA:C:N+1'].describe(), '', True)

    rep_mak.printReport()


