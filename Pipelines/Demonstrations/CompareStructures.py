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
recreate_csv,apply_dssp,recreate_html = 1,0,1

struc1dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
struc1 = ['7vos']
struc2dir = 'C:/Dev/Github/ProteinDataFiles/pdb_adj/'
struc2 = ['7vos']

csv_finala = "Csv/StructureCompare7vosA.csv"
csv_finalb = "Csv/StructureCompare7vosB.csv"
html_filename = 'Html/StructureCompare7vos.html'
geos = ['C:O','CA:C:N+1','C-1:N:CA:C','N:CA:C:N+1','C:N+1','N:CA:C']
geosX = ['CA:C:N+1','C:N+1']
geosY = ['C:O']
title = 'Two Structure Compare' # used at the header of the html report
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
import pandas as pd
###############################################################################################
if recreate_csv:
    print('### Load structures from BioPython #############')
    strucsa = bpm.loadPdbStructures(struc1,struc1dir,extension='ent',prefix='pdb',log=2)
    strucsb = bpm.loadPdbStructures(struc2, struc2dir, extension='ent', prefix='pdb', log=2)
    print('### Creating dataframe for correlations #############')
    geo_maka = dfm.GeometryMaker(strucsa, log=1)
    dataa = geo_maka.calculateGeometry(geos, log=1)
    geo_makb = dfm.GeometryMaker(strucsb, log=1)
    datab = geo_makb.calculateGeometry(geos, log=1)
    print('### Save dataframe ###############################')
    dataa.to_csv(csv_finala , index=False)
    print("Saved to", csv_finala)
    datab.to_csv(csv_finalb, index=False)
    print("Saved to", csv_finalb)

if apply_dssp:
    dataa = pd.read_csv(csv_finala + '_dssp.csv')
    datab = pd.read_csv(csv_finalb + '_dssp.csv')
else:
    dataa = pd.read_csv(csv_finala)
    datab = pd.read_csv(csv_finalb)

if recreate_html:
    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=2)
    # Some scatter plots are added to make better sense of the data, the colour scheme is to approcimate the AlphaFold probabilty
    rep_mak.addLineComment('Scatter Plots')

    for geox in geosX:
        for geoy in geosY:
            rep_mak.addPlot2d(dataa,'scatter',geox,geoy,hue='bfactor',title='Pdb Coordinates',palette='jet')
            rep_mak.addPlot2d(datab, 'scatter', geox, geoy, hue='bfactor', title='ED Maxima', palette='jet')
            if apply_dssp:
                rep_mak.addPlot2d(dataa, 'seaborn', geox, geoy, hue='ss', title='Pdb Coordinates', palette='tab10')
                rep_mak.addPlot2d(datab, 'seaborn', geox, geoy, hue='ss', title='ED Maxima', palette='tab10')
    rep_mak.printReport()


