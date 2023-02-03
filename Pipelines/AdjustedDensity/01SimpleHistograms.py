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
recreate_csv,recreate_html = 0,1

struc1dir = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
struc2dir = 'C:/Dev/Github/ProteinDataFiles/pdb_out/LapDen_ADJ/Density/'

structuresCsv = "Csv/pdb100.csv"

csv_finala = "Csv/StructureCompareA.csv"
csv_finalb = "Csv/StructureCompareB.csv"
html_filename = 'Html/StructureCompare.html'

titleA = "PDB Coordinates"
titleB = "Density Adjusted Coordinates"

geos = ['C:O','CA:C:N+1','C-1:N:CA:C','N:CA:C:N+1','C:N+1','N:CA:C',"CA:C","N:CA"]
geoHist = ["N:CA","CA:C",'C:O',"C:N+1"]
title = 'PDB vs Adjusted Coordinates Comparison'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm
from LeucipPy import HtmlReportMaker as hrm
import pandas as pd
###############################################################################################
if recreate_csv:
    print('### Load structures from BioPython #############')
    strucsdf = pd.read_csv(structuresCsv)
    strucs = strucsdf["PDB"].values


    strucsa = bpm.loadPdbStructures(strucs,struc1dir,extension='ent',prefix='pdb',log=2)
    strucsb = bpm.loadPdbStructures(strucs, struc2dir, extension='ent', prefix='pdb', log=2)
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

if recreate_html:

    dataa = pd.read_csv(csv_finala)
    datab = pd.read_csv(csv_finalb)

    # restrictions
    dataa = dataa.query("bfactor < 10")
    datab = datab.query("bfactor < 10")
    dataa = dataa.query("aa != 'SEP'")
    datab = datab.query("aa != 'SEP'")
    dataa = dataa.query("aa != 'STA'")
    datab = datab.query("aa != 'STA'")
    dataa = dataa.query("aa != 'IAS'")
    datab = datab.query("aa != 'IAS'")
    dataa = dataa.query("`aa+1` != 'SEP'")
    datab = datab.query("`aa+1` != 'SEP'")

    print(dataa.query("`CA:C` > 3")[["pdb_code","CA:C","infoCA:C"]])
    print(datab.query("`CA:C` > 3")[["pdb_code","CA:C","infoCA:C"]])

    print(dataa.query("`C:N+1` > 3")[["pdb_code", "C:N+1", "infoC:N+1"]])
    print(datab.query("`C:N+1` > 3")[["pdb_code", "C:N+1", "infoC:N+1"]])

    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=4)
    # Some scatter plots are added to make better sense of the data, the colour scheme is to approcimate the AlphaFold probabilty

    rep_mak.addLineComment(titleA)
    for geox in geoHist:
        rep_mak.addPlot1d(dataa,'histogram',geo_x=geox,title=geox, palette='steelblue',xrange=[1,1.6],bins=50)
    for geox in geoHist:
        rep_mak.addSeries(dataa[geox].describe(), '', True)
    rep_mak.addLineComment(titleB)
    for geox in geoHist:
        rep_mak.addPlot1d(datab,'histogram',geo_x=geox,title=geox, palette='crimson',xrange=[1,1.6],bins=50)
    for geox in geoHist:
        rep_mak.addSeries(datab[geox].describe(), '', True)

    rep_mak.printReport()


