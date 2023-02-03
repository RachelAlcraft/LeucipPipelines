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

csv_finala = "Csv/StructureCompareTauA.csv"
csv_finalb = "Csv/StructureCompareTauB.csv"
html_filename = 'Html/StructureCompareTau.html'

titleA = "PDB Coordinates"
titleB = "Density Adjusted Coordinates"

geos = ['C:O','CA:C:N+1','C-1:N:CA:C','N:CA:C:N+1','C:N+1','N:CA:C',"CA:C","N:CA","N:N+1","C-1:C"]
title = 'Bimodal Tau'
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
    dataa = dataa.query("aa != 'GLY'")
    datab = datab.query("aa != 'GLY'")
    dataa = dataa.query("aa != 'PRO'")
    datab = datab.query("aa != 'PRO'")

    print(dataa.query("`CA:C` > 3")[["pdb_code","CA:C","infoCA:C"]])
    print(datab.query("`CA:C` > 3")[["pdb_code","CA:C","infoCA:C"]])

    print(dataa.query("`N:N+1` > 4")[["pdb_code", "N:N+1", "infoN:N+1"]])
    print(datab.query("`N:N+1` > 4")[["pdb_code", "N:N+1", "infoN:N+1"]])

    dataa = dataa.query("`N:N+1` < 5")
    datab = datab.query("`N:N+1` < 5")


    print('### Creating html reports')
    rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=3)
    # Some scatter plots are added to make better sense of the data, the colour scheme is to approcimate the AlphaFold probabilty

    rep_mak.addLineComment(titleA)
    geox,geoy,hue = "N:CA:C","N:CA:C:N+1","C-1:N:CA:C"
    rep_mak.addPlot2d(dataa, 'scatter', geox, geoy, hue=hue, title=geox+"|"+geoy+"|"+hue, palette='jet')
    rep_mak.addPlot2d(dataa, 'probability', geox, geoy, hue=hue, title=geox+"|"+geoy+"|"+hue, palette="magma_r")
    geox, geoy, hue = "N:CA:C:N+1", "N:N+1", "N:CA:C"
    rep_mak.addPlot2d(dataa, 'scatter', geox, geoy, hue=hue, title=geox + "|" + geoy + "|" + hue, palette='jet')
    #rep_mak.addPlot2d(dataa, 'probability', geox, geoy, hue=hue, title=geox + "|" + geoy + "|" + hue, palette="inferno")

    geox, geoy, hue = "N:CA:C", "C-1:N:CA:C", "N:CA:C:N+1"
    rep_mak.addPlot2d(dataa, 'scatter', geox, geoy, hue=hue, title=geox+"|"+geoy+"|"+hue, palette='jet')
    rep_mak.addPlot2d(dataa, 'probability', geox, geoy, hue=hue, title=geox + "|" + geoy + "|" + hue, palette="magma_r")

    geox, geoy, hue = "C-1:N:CA:C", "C-1:C", "N:CA:C"
    rep_mak.addPlot2d(dataa, 'scatter', geox, geoy, hue=hue, title=geox + "|" + geoy + "|" + hue, palette='jet')
    #rep_mak.addPlot2d(dataa, 'probability', geox, geoy, hue=hue, title=geox + "|" + geoy + "|" + hue, palette="inferno")



    rep_mak.printReport()


