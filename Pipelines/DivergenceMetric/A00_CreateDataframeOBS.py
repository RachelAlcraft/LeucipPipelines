import pandas as pd
from Bio.PDB import PDBIO

from LeucipPy import HtmlReportMaker as hm
import A0Class_AtomCoordinates as coords
import Ext_Geometry as ext_geo
import Ext_PeptideBuilder as ext_pep
import A0_Globals as globals
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import GeometryMaker as dfm

def runMakeCsv(ID,PdbFile,PdbDir,AllGeos,aa_filter):
    print('LeucipPipeline: - Creating data -----------------------------------------------------------')
    pdb_df = pd.read_csv(PdbFile)
    pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
    pdbs = pdb_df['PDB'].values
    #pdbs = pdbs[:10]
    ### 2) Loop through and create bio python objects
    print('### Load structures from BioPython #############')
    strucs = bpm.loadPdbStructures(pdbs, PdbDir, extension='ent', prefix='pdb', log=2)

    print('### Creating dataframe for correlations #############')
    geo_mak = dfm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(AllGeos, log=1)

    print('### Save dataframe ###############################')
    print('filter on', aa_filter)
    if aa_filter != {}:
        print('filter on',aa_filter,len(data.index))
        data = geo_mak.filterDataFrame(data, inclusions={'aa':aa_filter})
        data = geo_mak.filterDataFrame(data, inclusions={'aa+1': aa_filter})
        data = geo_mak.filterDataFrame(data, inclusions={'aa-1': aa_filter})
        print('filterred', len(data.index))

    data = data.query('bfactor <= 10')
    data = data.query('occupancy == 1')

    print('Save to', "Csv/" + ID + "_01_Geometry.csv")
    data.to_csv("Csv/" + ID + "_01_Geometry.csv", index=False)
    return data

def makeHtmlReport(orig,id):
    geos = [] #for the display
    geos += ['N:CA', 'CA:C', 'C:O', 'N:C-1']  # ,'C:N+1']
    geos += ['C-1:N:CA', 'N:CA:C', 'CA-1:C-1:N', 'CA:C:O']
    geos += ['N-1:CA-1:C-1:N', 'N:CA:C:N+1', 'C-1:N:CA:C', 'N:CA-1:C-1:O-1']
    geos += ['N:O', 'O:C-1', 'N:N+1', 'O-1:N+1']
    geos += ['N:CA:C:N+1', 'CA:C:N+1:CA+1']

    print('### Creating report on synthetic structures #############')

    rep_mak = hm.HtmlReportMaker('Synthetic Structures Geometry', 'Html/AO_Observed_'+id+'.html', cols=4)
    rep_mak.addPlot2d(orig, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Observed Ramachandran')
    rep_mak.addPlot2d(orig, 'scatter', 'C-1:N:CA:C', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed Ramachandran-1')
    rep_mak.addPlot2d(orig, 'scatter', 'N:CA:C:N+1', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed psis')
    rep_mak.addPlot2d(orig, 'scatter', 'N:CA:C:N+1', 'N:O', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed psis')

    rep_mak.addLineComment('Bond lengths')
    for i in range(0,4):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='',overlay=False,alpha=0.8,density=True,bins=50)
    for i in range(0,4):# this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(),'Obs ' + geos[i],True)

    rep_mak.addLineComment('Bond angles')
    for i in range(4, 8):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='', overlay=False, alpha=0.8, density=True,bins=50)
    for i in range(4, 8):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(), 'Obs ' + geos[i], True)

    rep_mak.addLineComment('Dihedrals/Impropers')
    for i in range(8, 12):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='', overlay=False, alpha=0.8, density=True,bins=50)
    for i in range(8, 12):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(), 'Obs ' + geos[i], True)

    rep_mak.addLineComment('Some other paramaters')
    for i in range(12, 16):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='', overlay=False, alpha=0.8, density=True,bins=50)
    for i in range(12, 16):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(), 'Obs ' + geos[i], True)

    rep_mak.printReport()

################ MAIN CODE RUN HERE ##################
#IDS = ['High_GLY','Redo_GLY']
#for ID in IDS:
import A0_Globals as globals
from datetime import datetime
import sys
sys.path.append('../1Library')
import Helpers as hlp

PdbFile = 'C:/Dev/Github/LeucipPipelines/Pipelines/0ConfigData/Pdbs_Under1_good_redo.csv'
PdbDirHigh = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
PdbDirRedo = "C:/Dev/Github/ProteinDataFiles/pdb_data_redo/"
ID_high = 'PW_High'
ID_redo = 'PW_Redo'

runs = [1]
for run in runs:
    start = datetime.now()
    geos = globals.getGeos(False,False)
    geosGLY = globals.getGeos(True,False)
    for run_for in runs:
        print(run_for)
        if 1 == run_for:
            data = runMakeCsv(ID_high,PdbFile,PdbDirHigh,geos,aa_filter=globals.aa_NO_gly_pro)
            makeHtmlReport(data,ID_high)
        if 2 == run_for:
            data = runMakeCsv(ID_high+'_GLY', PdbFile, PdbDirHigh, geosGLY,aa_filter=globals.aa_GLY)
            makeHtmlReport(data,ID_high+'_GLY')
        if 3 == run_for:
            data = runMakeCsv(ID_redo, PdbFile, PdbDirRedo, geos,aa_filter=globals.aa_NO_gly_pro)
            makeHtmlReport(data,ID_redo)
        if 4 == run_for:
            data = runMakeCsv(ID_redo+'_GLY', PdbFile, PdbDirRedo, geosGLY,aa_filter=globals.aa_GLY)
            makeHtmlReport(data,ID_redo+'_GLY')
    end = datetime.now()
    hlp.printTime(start,end)
