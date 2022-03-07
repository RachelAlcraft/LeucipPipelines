import pandas as pd
from Bio.PDB import PDBIO

from LeucipPy import GeometryMaker as gm
from LeucipPy import HtmlReportMaker as hm
import A0Class_AtomCoordinates as coords
import Ext_Geometry as ext_geo
import Ext_PeptideBuilder as ext_pep
import A0_Globals as glob


def makePdbSynthetic(data,numStrucs,lenChain,pair_rama):
    print('LeucipPipeline: - Creating synthetic data -----------------------------------------------------------')
    strucs = []
    last_psi = None
    last_improper = None
    for i in range(numStrucs):
        geo = ext_geo.geometry('G')
        structure = ext_pep.initialize_res(geo)
        structure.header = {}
        structure.header['resolution'] = 1
        structure.id = 'X' + str(i)
        for j in range(lenChain):
            param_gen = coords.AtomCoordinates(data, log=1)
            if pair_rama == 1:
                geom, psi_next,improper_next = param_gen.generateRamaGeo()
            elif pair_rama == 2:
                geom,psi_next,improper_next = param_gen.generateSampleGeo()
            else:
                geom,psi_next,improper_next = param_gen.generateAllRandomGeo()

            if last_psi != None:
                geom.psi_im1 = last_psi
                #geom.N_CA_C_O_diangle = last_improper which one?
            structure = ext_pep.add_residue(structure, geom)
            last_psi = psi_next
            last_improper = improper_next
        strucs.append(structure)
    return strucs

def makeGeometryCsv(strucs,pair_rama):
    print('Make geo csv')
    geos = glob.getGeos(True)
    print(geos)
    print('### Creating dataframe for correlations #############')
    geo_mak = gm.GeometryMaker(strucs, log=1)
    data = geo_mak.calculateGeometry(geos, log=1)
    data = glob.trimData(data,15,geos)

    csv_file = "Csv/" + ID + "_IDEAL_UNPAIRED_01_Geometry.csv"
    if pair_rama == 1:
        csv_file = "Csv/" + ID + "_IDEAL_PAIRED_01_Geometry.csv"
    elif pair_rama == 2:
        csv_file = "Csv/" + ID + "_IDEAL_SAMPLED_01_Geometry.csv"

    print('Save to', csv_file)
    data.to_csv(csv_file, index=False)
    print(data)
    return data

def makeHtmlReport(ID,data,orig,pair_rama):
    geos = [] #for the display
    geos += ['N:CA', 'CA:C', 'C:O', 'N:C-1']  # ,'C:N+1']
    geos += ['C-1:N:CA', 'N:CA:C', 'CA-1:C-1:N', 'CA:C:O']
    geos += ['N-1:CA-1:C-1:N', 'N:CA:C:N+1', 'C-1:N:CA:C', 'N:CA-1:C-1:O-1']
    geos += ['N:O', 'O:C-1', 'N:N+1', 'O-1:N+1']
    geos += ['N:CA:C:N+1', 'CA:C:N+1:CA+1']

    # some need to be abs because otherwise they are hard to visualise
    data['N:CA-1:C-1:O-1'] = abs(data['N:CA-1:C-1:O-1'])
    orig['N:CA-1:C-1:O-1'] = abs(orig['N:CA-1:C-1:O-1'])

    print('### Creating report on synthetic structures #############')

    rep_mak = hm.HtmlReportMaker('Synthetic Structures Geometry', 'Html/AO_' + ID + '_SyntheticCoords_'+str(pair_rama)+'.html', cols=4)
    rep_mak.addPlot2d(orig, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Observed Ramachandran')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Synthetic Ramachandran')
    rep_mak.addPlot2d(orig, 'scatter', 'C-1:N:CA:C', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed Ramachandran-1')
    rep_mak.addPlot2d(data, 'scatter', 'C-1:N:CA:C', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral', title='Synthetic Ramachandran-1')

    rep_mak.addPlot2d(orig, 'scatter', 'N:CA:C:N+1', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed psis')
    rep_mak.addPlot2d(data, 'scatter', 'N:CA:C:N+1', 'N-1:CA-1:C-1:N', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Synthetic psis')
    rep_mak.addPlot2d(orig, 'scatter', 'N:CA:C:N+1', 'N:O', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Observed psis')
    rep_mak.addPlot2d(data, 'scatter', 'N:CA:C:N+1', 'N:O', hue='N-1:CA-1:C-1:N', palette='Spectral',title='Synthetic psis')

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
IDS = ['PW_High_GLY','PW_High']
for ID in IDS:
    print('Load from', "Csv/" + ID + "_01_Geometry.csv")
    data = pd.read_csv("Csv/" + ID + "_01_Geometry.csv")
    for pair_rama in [0, 1]:
        strucs = makePdbSynthetic(data, 5000,6,pair_rama)
        io = PDBIO()
        io.set_structure(strucs[0])
        io.save('Csv/syn' + str(pair_rama) + '.pdb')
        csv = makeGeometryCsv(strucs,pair_rama)
        makeHtmlReport(ID,csv,data,pair_rama)
