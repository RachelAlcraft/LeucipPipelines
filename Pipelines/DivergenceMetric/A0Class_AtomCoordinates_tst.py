import pandas as pd
from Bio.PDB import PDBIO

from LeucipPy import GeometryMaker as gm
from LeucipPy import HtmlReportMaker as hm
import A0Class_AtomCoordinates as coords
import Ext_Geometry as ext_geo
import Ext_PeptideBuilder as ext_pep
geos = []
geos += ['N:CA','CA:C','C:O','N:C-1']#,'C:N+1']
geos += ['C-1:N:CA','N:CA:C','CA-1:C-1:N','CA:C:O']
geos += ['N-1:CA-1:C-1:N', 'C-1:N:CA:C', 'CA-1:C-1:N:CA','N:CA:C:O']
geos += ['N:O', 'O:C-1', 'N:N+1','O-1:N+1']
geos += ['N:CA:C:N+1','CA:C:N+1:CA+1']

def makePdbSynthetic(data,numStrucs,lenChain):

    print('LeucipPipeline: - Creating synthetic data -----------------------------------------------------------')
    strucs = []
    last_psi = None
    for i in range(numStrucs):
        geo = ext_geo.geometry('G')
        structure = ext_pep.initialize_res(geo)
        structure.header = {}
        structure.header['resolution'] = 1
        structure.id = 'X' + str(i)
        for j in range(lenChain):
            param_gen = coords.AtomCoordinates(data, log=1)
            #geom,psi_next = param_gen.generateAllRandomGeo()
            geom,psi_next = param_gen.generateRamaGeo()
            #geom,psi_next = param_gen.generateSampleGeo()
            if last_psi != None:
                #print('diff1', geom.psi_im1, last_psi)
                geom.psi_im1 = last_psi
            structure = ext_pep.add_residue(structure, geom)
            print('diff2', geom.psi_im1, last_psi)
            last_psi = psi_next
        strucs.append(structure)
    return strucs

def makeGeometryCsv(strucs):
    print('### Creating dataframe for correlations #############')
    geo_mak = gm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(geos, log=1)
    print('Save to', "Csv/PW_TSTL_01_Geometry.csv")
    data.to_csv("Csv/PW_TST_01_Geometry.csv", index=False)
    return data

def makeHtmlReport(data,orig):
    print('### Creating report on synthetic structures #############')
    #make 2 of them abs val
    #data['CA-1:C-1:N'] = abs(data['CA:C:N+1:CA+1'])
    #orig['CA:C:N+1'] = abs(orig['CA:C:N+1:CA+1'])
    #data['CA:C:O:N+1'] = abs(data['CA:C:O:N+1'])
    #orig['CA:C:O:N+1'] = abs(orig['CA:C:O:N+1'])

    rep_mak = hm.HtmlReportMaker('Synthetic Structures Geometry', 'Html/AO_SyntheticRandomCoords.html', cols=4)
    rep_mak.addPlot2d(orig, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Observed Ramachandran')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Synthetic Ramachandran')
    rep_mak.addPlot2d(orig, 'scatter', 'C-1:N:CA:C', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed Ramachandran-1')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Synthetic Ramachandran-1')

    rep_mak.addPlot2d(orig, 'scatter', 'N:CA:C:N+1', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed psis')
    rep_mak.addPlot2d(data, 'scatter', 'N:CA:C:N+1', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Synthetic psis')

    #rep_mak.addPlot2d(orig, 'scatter', 'N:CA', 'CA:C', hue='C:O', palette='Spectral', title='Observed')
    #rep_mak.addPlot2d(data, 'scatter', 'N:CA', 'CA:C', hue='C:O', palette='Spectral',title='Synthetic')

    rep_mak.addLineComment('Comparing bond lengths that were used in synthetic creation')
    for i in range(0,4):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='',overlay=True,alpha=0.8,density=True,bins=50)
        rep_mak.addPlot1d(data, 'histogram', geos[i], hue='aa', title='Syn=blue, Obs=red', alpha=0.5, density=True,palette='navy',bins=50)
    for i in range(0,4):# this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(),'Obs ' + geos[i],True)
    for i in range(0, 4):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(data[geos[i]].describe(), 'Syn ' + geos[i], True)

    rep_mak.addLineComment('Comparing bond angles that were used in synthetic creation')
    for i in range(4, 8):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='', overlay=True, alpha=0.8, density=True,bins=50)
        rep_mak.addPlot1d(data, 'histogram', geos[i], hue='aa', title='Syn=blue, Obs=red', alpha=0.5, density=True, palette='navy',bins=50)
    for i in range(4, 8):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(), 'Obs ' + geos[i], True)
    for i in range(4, 8):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(data[geos[i]].describe(), 'Syn ' + geos[i], True)

    rep_mak.addLineComment('Comparing bond dihedrals that were used in synthetic creation')
    for i in range(8, 12):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='', overlay=True, alpha=0.8, density=True,bins=50)
        rep_mak.addPlot1d(data, 'histogram', geos[i], hue='aa', title='Syn=blue, Obs=red', alpha=0.5, density=True, palette='navy',bins=50)
    for i in range(8, 12):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(), 'Obs ' + geos[i], True)
    for i in range(8, 12):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(data[geos[i]].describe(), 'Syn ' + geos[i], True)

    rep_mak.addLineComment('Comparing some paramaters not used in the creation')
    for i in range(12, 16):
        rep_mak.addPlot1d(orig, 'histogram', geos[i], hue='aa', title='', overlay=True, alpha=0.8, density=True,bins=50)
        rep_mak.addPlot1d(data, 'histogram', geos[i], hue='aa', title='Syn=blue, Obs=red', alpha=0.5, density=True, palette='navy',bins=50)
    for i in range(12, 16):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(orig[geos[i]].describe(), 'Obs ' + geos[i], True)
    for i in range(12, 16):  # this adds the summary statistics below each of the above histgrams
        rep_mak.addSeries(data[geos[i]].describe(), 'Syn ' + geos[i], True)

    rep_mak.printReport()

################ MAIN CODE RUN HERE ##################
ID = 'High_GLY'
print('Load from', "Csv/PW_" + ID + "_01_Geometry.csv")
data = pd.read_csv("Csv/PW_" + ID + "_01_Geometry.csv")
strucs = makePdbSynthetic(data, 1000,3)
io = PDBIO()
io.set_structure(strucs[0])
io.save('Csv/syn.pdb')
csv = makeGeometryCsv(strucs)
makeHtmlReport(csv,data)
