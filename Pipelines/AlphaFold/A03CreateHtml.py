import sys

from LeucipPy import GeoHTML as ghm
from LeucipPy import GeoDataFrame as gdf
import pandas as pd
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Helpers as help

def compareProteomes(dfs):
    html_name = 'Html/AlphaFold_CompareProteomes.html'
    title = "Compare Proteomes Geometry from AlphaFold"

    georep = ghm.GeoHTML(title, html_name, cols=4)
    for name,data,data_corr in dfs:
        print(name + ' Proteome ascending-------------------')
        #data['Probability'] = data['bfactor']
        data = data.sort_values(by='Probability', ascending=True)
        georep.addLineComment(name + ' Proteome')
        georep.addPlot2d(data=data, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title=name + ' Proteome', crange=[50, 90])
        georep.addPlot2d(data=data, plottype='scatter', geo_x='N:CA', geo_y='CA:C', hue='Probability', alpha=0.6, palette='jet_r', title=name + ' Proteome', crange=[50, 90])
        georep.addPlot2d(data=data, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title=name + ' Proteome', crange=[50, 90])
        georep.addPlot2d(data=data, plottype='scatter', geo_x='N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title=name + ' Proteome', crange=[50, 90])

    for name,data,data_corr in dfs:
        print(name + ' Proteome descending-------------------')
        #data['Probability'] = data['bfactor']
        #data = data.sort_values(by='Probability', ascending=True)
        georep.addLineComment(name + ' Proteome >=90% probable')
        data = data.query('Probability >= 80')
        georep.addPlot2d(data=data, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='Probability', alpha=0.6, palette='viridis', title=name + ' Proteome', crange=[90, 100])
        georep.addPlot2d(data=data, plottype='scatter', geo_x='N:CA', geo_y='CA:C', hue='Probability', alpha=0.6, palette='viridis', title=name + ' Proteome', crange=[90, 100])
        georep.addPlot2d(data=data, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability',alpha=0.6, palette='viridis', title=name + ' Proteome', crange=[90, 100])
        georep.addPlot2d(data=data, plottype='scatter', geo_x='N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='viridis', title=name + ' Proteome', crange=[90, 100])

    for name,data,data_corr in dfs:
        print(name + ' Proteome correlations-------------------')
        #data_corr['Probability'] = data_corr['bfactor']
        data_corr = data_corr.sort_values(by='Probability', ascending=True)
        georep.addLineComment(name + ' Proteome correlations >=90% probable')
        data_corr = data_corr.query('Probability >= 90')
        georep.addPlot2d(data=data_corr, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='N:CA:C', alpha=0.8,palette='rainbow', title='AlphaFold')
        georep.addPlot2d(data=data_corr, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='N:CA:C', alpha=0.8,palette='rainbow', title='AlphaFold')
        georep.addPlot2d(data=data_corr, plottype='scatter', geo_x='N:O', geo_y='CB:O', hue='Probability', alpha=0.6,palette='viridis', title='AlphaFold', crange=[90, 100])
        georep.addPlot2d(data=data_corr, plottype='scatter', geo_x='CA-2:CA-1:CA', geo_y='CA:CA+1:CA+2',hue='Probability', alpha=0.6, palette='viridis', title='AlphaFold', crange=[90, 100])


    georep.changeColNumber(1)
    geos = [['N:CA', 1,2], ['CA:C', 1,2], ['C:O', 1,2], ['C:N+1', 1,2], ['N:CA:C',90, 130]]
    for name,data,data_corr in dfs:
        print(name + ' Proteome extremes -------------------')
        georep.addLineComment('Extremes for proteome >= 90% ' +name)
        #data['Probability'] = data['bfactor']
        data['Link'] = data.apply(lambda row: help.applyALPHALINK(row['pdb_code']), axis=1)
        data = data.sort_values(by='Probability', ascending=False)
        for geo,lower,upper in geos:
            data = data.query('Probability >= 90')
            dataExtUp = data.query('`'+geo+'` >=' + str(upper))
            dataExtLow = data.query('`' + geo + '` <=' + str(lower))
            if len(dataExtUp.index) > 0:
                dataExtUp = dataExtUp[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability','Link']]
                georep.addDataFrame(dataExtUp)
            if len(dataExtLow.index) > 0:
                dataExtLow = dataExtLow[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability','Link']]
                georep.addDataFrame(dataExtLow)

    georep.printReport()
    print("Saved report to",html_name)

def createHtml(name,data,datacorr,geoChecksA,geoChecksB):
    #data['Probability'] = data['bfactor']
    #datacorr['Probability'] = datacorr['bfactor']
    data['Link'] = data.apply(lambda row: help.applyALPHALINK(row['pdb_code']),axis=1)
    datacorr['Link'] = datacorr.apply(lambda row: help.applyALPHALINK(row['pdb_code']), axis=1)

    print('Creating sorted dataframes...')
    print('...probability ascending')
    probasc = data.sort_values(by='Probability', ascending=True)
    print('...probability descending')
    probdesc = data.sort_values(by='Probability', ascending=False)
    print('...correlation probability ascending')
    probcorrasc = datacorr.sort_values(by='Probability', ascending=True)
    print('...correlation probability descending')
    probcorrdesc = datacorr.sort_values(by='Probability', ascending=False)
    print('...correlation tau ascending')
    tauasc = datacorr.sort_values(by='N:CA:C', ascending=True)

    dataTauProb = tauasc.query('Probability >= 90')
    dataGeoProb = probasc.query('Probability >= 90')

    html_name = 'Html/AlphaFold_' + name + '.html'
    title = name + ": Explore Geometry from AlphaFold"
    print(name, ' Creating HTML files #############################')
    georep = ghm.GeoHTML(title, html_name, cols=3)
    georep.changeColNumber(5)
    print(' ...  Line 1 sort  #############################')
    georep.addLineComment(name + ' Proteome, most probable at the top')
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='N:CA', geo_y='CA:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='CA:C:N+1', alpha=0.6, palette='rainbow', title='All probability')
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])


    print(' ...  Line 2 reverse sort  #############################')
    georep.addLineComment(name + ' Proteome, least probable at the top')
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='N:CA', geo_y='CA:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=dataGeoProb, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='CA:C:N+1', alpha=0.6, palette='rainbow', title='>= 90')
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability',alpha=0.6, palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])

    georep.changeColNumber(4)
    print(' ...  Line 3 correlations  #############################')
    georep.addLineComment(name + ' Proteome Geometric Correlations, most probable at the top')
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='N:O', geo_y='CB:O', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='CA-2:CA-1:CA', geo_y='CA:CA+1:CA+2', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])

    print(' ...  Line 4 correlations  #############################')
    datacorr = datacorr.sort_values(by='Probability', ascending=False)
    georep.addLineComment(name + ' Proteome Geometric Correlations, least probable at the top')
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='N:O', geo_y='CB:O', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='CA-2:CA-1:CA', geo_y='CA:CA+1:CA+2', hue='Probability',alpha=0.6, palette='jet_r', title='', crange=[50, 90])

    print(' ...  Line 4 correlations  #############################')
    georep.addLineComment(name + ' Proteome Geometric Correlations, hue on tau')
    georep.addPlot2d(data=tauasc, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='N:CA:C', alpha=0.6,palette='rainbow', title='All probability',xrange=[-180,180])
    georep.addPlot2d(data=tauasc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='N:CA:C', alpha=0.6,palette='rainbow', title='All probability',xrange=[-180,180])
    georep.addPlot2d(data=dataTauProb, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='N:CA:C', alpha=0.6,palette='rainbow', title='>= 90%',xrange=[-180,180])
    georep.addPlot2d(data=dataTauProb, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='N:CA:C', alpha=0.6,palette='rainbow', title='>= 90%',xrange=[-180,180])


    #filter extremes
    georep.changeColNumber(1)
    georep.addLineComment('Extreme backbone values >= 80%')
    for geo, lower, upper in geoChecksA:
        print(geo,'extremes  #############################')
        dataExt = data.query('`'+geo+'` >=' + str(upper))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Upper extreme values,'+ geo)
            dataExt = dataExt[['pdb_code',geo,'chain','rid','aa','Probability','info'+geo,'Link']]
            georep.addDataFrame(dataExt)
        print(geo, 'extremes  #############################')
        dataExt = data.query('`' + geo + '` <=' + str(lower))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Lower extreme values,'+ geo)
            dataExt = dataExt[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability', 'info' + geo,'Link']]
            georep.addDataFrame(dataExt)

    for geo, lower, upper in geoChecksB:
        print(geo, 'extremes  #############################')
        dataExt = datacorr.query('`' + geo + '` >=' + str(upper))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Upper extreme values,'+ geo)
            dataExt = dataExt[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability', 'info' + geo,'Link']]
            georep.addDataFrame(dataExt)
        print(geo, 'extremes  #############################')
        dataExt = datacorr.query('`' + geo + '` <=' + str(lower))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Lower extreme values,'+ geo)
            dataExt = dataExt[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability', 'info' + geo,'Link']]
            georep.addDataFrame(dataExt)

    georep.printReport()
    print("Saved report to",html_name)

def createHtmlFromGeos(name,ext,data,geosone, geostwo):
    data['Link'] = data.apply(lambda row: help.applyALPHALINK(row['pdb_code']), axis=1)
    html_name = 'Html/AlphaFold_' + name + "_" + ext + '.html'
    title = name + ": Explore Geometry from AlphaFold"
    print(name, ' Creating HTML files #############################')

    georep = ghm.GeoHTML(title, html_name, cols=3)
    georep.changeColNumber(len(geosone))
    print(' ...  Line 1 histograms  #############################')
    georep.addLineComment(name + ' Proteome, histograms')
    for geo in geosone:
        georep.addPlot1d(data,'histogram',geo)

    georep.changeColNumber(len(geostwo))
    georep.addLineComment(name + ' Proteome, scatters')
    for geox,geoy in geostwo:
        print(geox,geoy,'correlations...')
        data = data.sort_values(by='Probability', ascending=True)
        georep.addPlot2d(data,'scatter',geox,geoy, hue='Probability', alpha=0.6,palette='jet_r', title=name + 'AlphaFold', crange=[50, 90])

    #filter extremes
    georep.changeColNumber(1)
    for geo in geosone:
        print(' ...',geo,'extremes...')
        data = data.sort_values(by=geo, ascending=False)
        dataMax = data[:10]
        dataMin = data[-10:]
        georep.addLineComment('Extreme ' +  geo + ' values')
        dataMax = dataMax[['pdb_code',geo,'chain','rid','aa','Probability','info'+geo]]
        dataMin = dataMin[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability','info'+geo]]
        georep.addDataFrame(dataMin)
        georep.addDataFrame(dataMax)



    georep.printReport()
    print("Saved report to",html_name)


def createHtml(name,data,datacorr,geoChecksA,geoChecksB):
    #data['Probability'] = data['bfactor']
    #datacorr['Probability'] = datacorr['bfactor']
    data['Link'] = data.apply(lambda row: help.applyALPHALINK(row['pdb_code']),axis=1)
    datacorr['Link'] = datacorr.apply(lambda row: help.applyALPHALINK(row['pdb_code']), axis=1)

    print('Creating sorted dataframes...')
    print('...probability ascending')
    probasc = data.sort_values(by='Probability', ascending=True)
    print('...probability descending')
    probdesc = data.sort_values(by='Probability', ascending=False)
    print('...correlation probability ascending')
    probcorrasc = datacorr.sort_values(by='Probability', ascending=True)
    print('...correlation probability descending')
    probcorrdesc = datacorr.sort_values(by='Probability', ascending=False)
    print('...correlation tau ascending')
    tauasc = datacorr.sort_values(by='N:CA:C', ascending=True)

    dataTauProb = tauasc.query('Probability >= 90')
    dataGeoProb = probasc.query('Probability >= 90')

    html_name = 'Html/AlphaFold_' + name + '.html'
    title = name + ": Explore Geometry from AlphaFold"
    print(name, ' Creating HTML files #############################')
    georep = ghm.GeoHTML(title, html_name, cols=3)
    georep.changeColNumber(5)
    print(' ...  Line 1 sort  #############################')
    georep.addLineComment(name + ' Proteome, most probable at the top')
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='N:CA', geo_y='CA:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='CA:C:N+1', alpha=0.6, palette='rainbow', title='All probability')
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probasc, plottype='scatter', geo_x='N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6, palette='jet_r', title='', crange=[50, 90])


    print(' ...  Line 2 reverse sort  #############################')
    georep.addLineComment(name + ' Proteome, least probable at the top')
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='N:CA', geo_y='CA:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=dataGeoProb, plottype='scatter', geo_x='C:O', geo_y='C:N+1', hue='CA:C:N+1', alpha=0.6, palette='rainbow', title='>= 90')
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='N:CA:C:N+1', hue='Probability',alpha=0.6, palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probdesc, plottype='scatter', geo_x='N:CA:C', geo_y='N:CA:C:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])

    georep.changeColNumber(4)
    print(' ...  Line 3 correlations  #############################')
    georep.addLineComment(name + ' Proteome Geometric Correlations, most probable at the top')
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='N:O', geo_y='CB:O', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrasc, plottype='scatter', geo_x='CA-2:CA-1:CA', geo_y='CA:CA+1:CA+2', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])

    print(' ...  Line 4 correlations  #############################')
    datacorr = datacorr.sort_values(by='Probability', ascending=False)
    georep.addLineComment(name + ' Proteome Geometric Correlations, least probable at the top')
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='N:O', geo_y='CB:O', hue='Probability', alpha=0.6,palette='jet_r', title='', crange=[50, 90])
    georep.addPlot2d(data=probcorrdesc, plottype='scatter', geo_x='CA-2:CA-1:CA', geo_y='CA:CA+1:CA+2', hue='Probability',alpha=0.6, palette='jet_r', title='', crange=[50, 90])

    print(' ...  Line 4 correlations  #############################')
    georep.addLineComment(name + ' Proteome Geometric Correlations, hue on tau')
    georep.addPlot2d(data=tauasc, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='N:CA:C', alpha=0.6,palette='rainbow', title='All probability',xrange=[-180,180])
    georep.addPlot2d(data=tauasc, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='N:CA:C', alpha=0.6,palette='rainbow', title='All probability',xrange=[-180,180])
    georep.addPlot2d(data=dataTauProb, plottype='scatter', geo_x='N:CA:C:N+1', geo_y='N:N+1', hue='N:CA:C', alpha=0.6,palette='rainbow', title='>= 90%',xrange=[-180,180])
    georep.addPlot2d(data=dataTauProb, plottype='scatter', geo_x='C-1:N:CA:C', geo_y='C-1:C', hue='N:CA:C', alpha=0.6,palette='rainbow', title='>= 90%',xrange=[-180,180])


    #filter extremes
    georep.changeColNumber(1)
    georep.addLineComment('Extreme backbone values >= 80%')
    for geo, lower, upper in geoChecksA:
        print(geo,'extremes  #############################')
        dataExt = data.query('`'+geo+'` >=' + str(upper))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Upper extreme values,'+ geo)
            dataExt = dataExt[['pdb_code',geo,'chain','rid','aa','Probability','info'+geo,'Link']]
            georep.addDataFrame(dataExt)
        print(geo, 'extremes  #############################')
        dataExt = data.query('`' + geo + '` <=' + str(lower))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Lower extreme values,'+ geo)
            dataExt = dataExt[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability', 'info' + geo,'Link']]
            georep.addDataFrame(dataExt)

    for geo, lower, upper in geoChecksB:
        print(geo, 'extremes  #############################')
        dataExt = datacorr.query('`' + geo + '` >=' + str(upper))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Upper extreme values,'+ geo)
            dataExt = dataExt[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability', 'info' + geo,'Link']]
            georep.addDataFrame(dataExt)
        print(geo, 'extremes  #############################')
        dataExt = datacorr.query('`' + geo + '` <=' + str(lower))
        dataExt = dataExt.query('Probability >= 80')
        if len(dataExt.index) > 0:
            georep.addLineComment('Lower extreme values,'+ geo)
            dataExt = dataExt[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability', 'info' + geo,'Link']]
            georep.addDataFrame(dataExt)

    georep.printReport()
    print("Saved report to",html_name)

def createHtmlFromGeos(name,ext,data,geosone, geostwo):
    data['Link'] = data.apply(lambda row: help.applyALPHALINK(row['pdb_code']), axis=1)
    html_name = 'Html/AlphaFold_' + name + "_" + ext + '.html'
    title = name + ": Explore Geometry from AlphaFold"
    print(name, ' Creating HTML files #############################')

    georep = ghm.GeoHTML(title, html_name, cols=3)
    georep.changeColNumber(len(geosone))
    print(' ...  Line 1 histograms  #############################')
    georep.addLineComment(name + ' Proteome, histograms')
    for geo in geosone:
        georep.addPlot1d(data,'histogram',geo)

    georep.changeColNumber(len(geostwo))
    georep.addLineComment(name + ' Proteome, scatters')
    for geox,geoy in geostwo:
        print(geox,geoy,'correlations...')
        data = data.sort_values(by='Probability', ascending=True)
        georep.addPlot2d(data,'scatter',geox,geoy, hue='Probability', alpha=0.6,palette='jet_r', title=name + 'AlphaFold', crange=[50, 90])

    #filter extremes
    georep.changeColNumber(1)
    for geo in geosone:
        print(' ...',geo,'extremes...')
        data = data.sort_values(by=geo, ascending=False)
        dataMax = data[:10]
        dataMin = data[-10:]
        georep.addLineComment('Extreme ' +  geo + ' values')
        dataMax = dataMax[['pdb_code',geo,'chain','rid','aa','Probability','info'+geo]]
        dataMin = dataMin[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability','info'+geo]]
        georep.addDataFrame(dataMin)
        georep.addDataFrame(dataMax)



    georep.printReport()
    print("Saved report to",html_name)


def createProbability(name, ext,df_geometry, df_correlation,scattersA, scattersB,allone=False):
    html_name = 'Html/AlphaFold_' + name + "_" + ext + '.html'
    title = name + ": Explore Geometry from AlphaFold"

    df_geometry['Probability'] = df_geometry['bfactor']
    df_correlation['Probability'] = df_correlation['bfactor']

    if not allone:
        print(name, ' Filtering DFs A #############################')
        df_geo_90 = df_geometry.query('Probability >= 90')
        df_geo_7090 = df_geometry.query('Probability < 90')
        df_geo_7090 = df_geo_7090.query('Probability >= 70')
        df_geo_5070 = df_geometry.query('Probability < 70')
        df_geo_5070 = df_geo_5070.query('Probability >= 50')
        df_geo_50 = df_geometry.query('Probability < 50')

        print(name, ' Filtering DFs B #############################')
        df_corr_90 = df_correlation.query('Probability >= 90')
        df_corr_7090 = df_correlation.query('Probability < 90')
        df_corr_7090 = df_corr_7090.query('Probability >= 70')
        df_corr_5070 = df_correlation.query('Probability < 70')
        df_corr_5070 = df_corr_5070.query('Probability >= 50')
        df_corr_50 = df_correlation.query('Probability < 50')
    else:
        print(name, ' Filtering DFs A #############################')
        df_geo_90 = df_geometry
        df_geo_7090 = df_geometry
        df_geo_7090 = df_geo_7090
        df_geo_5070 = df_geometry
        df_geo_5070 = df_geo_5070
        df_geo_50 = df_geometry

        print(name, ' Filtering DFs B #############################')
        df_corr_90 = df_correlation
        df_corr_7090 = df_correlation
        df_corr_7090 = df_corr_7090
        df_corr_5070 = df_correlation
        df_corr_5070 = df_corr_5070
        df_corr_50 = df_correlation

    print(name, ' Creating HTML files #############################')

    georep = ghm.GeoHTML(title, html_name, cols=4)
    print(' ...  Scatters A  #############################')
    georep.addLineComment(name + ' Proteome, probability plots')
    for geoX,geoY in scattersA:
        print('Creating >90',geoX,geoY)
        georep.addPlot2d(df_geo_90,'probability',geoX,geoY, hue='Probability', alpha=1,palette='cubehelix_r', title=name + 'AlphaFold >=90%')
        print('Creating 70-90', geoX, geoY)
        georep.addPlot2d(df_geo_7090, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r',title=name + 'AlphaFold 70-90%')
        print('Creating 50-70', geoX, geoY)
        georep.addPlot2d(df_geo_5070, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r',title=name + 'AlphaFold 50-70%')
        print('Creating <50', geoX, geoY)
        georep.addPlot2d(df_geo_50, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r',title=name + 'AlphaFold <50%')

    print(' ...  Scatters B  #############################')
    georep.addLineComment(name + ' Proteome, probability plots')
    for geoX, geoY in scattersB:
        print('Creating >90', geoX, geoY)
        if len(df_corr_90.index) > 0:
            georep.addPlot2d(df_corr_90, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r', title=name + 'AlphaFold >=90%')
        else:
            georep.addBoxComment('Empty for >90')
        print('Creating 70-90', geoX, geoY)
        if len(df_corr_7090.index) > 0:
            georep.addPlot2d(df_corr_7090, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r', title=name + 'AlphaFold 70-90%')
        else:
            georep.addBoxComment('Empty for 70-90')
        print('Creating 50-70', geoX, geoY)
        if len(df_corr_5070.index) > 0:
            georep.addPlot2d(df_corr_5070, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r', title=name + 'AlphaFold 50-70%')
        else:
            georep.addBoxComment('Empty for 50-70')
        print('Creating <50', geoX, geoY)
        if len(df_corr_50.index) > 0:
            georep.addPlot2d(df_corr_50, 'probability', geoX, geoY, hue='Probability', alpha=1, palette='cubehelix_r', title=name + 'AlphaFold <50%')
        else:
            georep.addBoxComment('Empty for <50')


    georep.printReport()
    print("Saved report to",html_name)

