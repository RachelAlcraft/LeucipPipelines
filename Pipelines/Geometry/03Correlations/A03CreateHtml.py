from LeucipPy import GeoHTML as ghm
from LeucipPy import GeoDataFrame as gdf
import pandas as pd


def createHtmlFromGeos(name,ext,data,geos):
    html_name = 'Html/AlphaFold_' + name + "_" + ext + '.html'
    title = name + ": Explore Geometry from AlphaFold"
    print(name, ' Creating HTML files #############################')
    data['Probability'] = data['bfactor']

    georep = ghm.GeoHTML(title, html_name, cols=3)
    georep.changeColNumber(len(geos))
    print(' ...  Line 1 histograms  #############################')
    georep.addLineComment(name + ' Proteome, histograms')
    for geo in geos:
        georep.addPlot1d(data,'histogram',geo)

    #filter extremes
    georep.changeColNumber(1)
    for geo in geos:
        print(' ...',geo,'extremes...')
        data = data.sort_values(by=geo, ascending=False)
        dataMax = data[:10]
        dataMin = data[-10:]
        georep.addLineComment('Extreme ' +  geo + ' values')
        dataMax = dataMax[['pdb_code',geo,'chain','rid','aa','Probability']]
        dataMin = dataMin[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability']]
        georep.addDataFrame(dataMin)
        georep.addDataFrame(dataMax)


    georep.printReport()
    print("Saved report to",html_name)

def createHtmlFromPairs(name,ext,data,pairs,geos):
    html_name = 'Html/' + name + "_" + ext + '_Scatters.html'
    title = name + ": Explore Geometry from " + name
    print(name, ' Creating HTML files #############################')
    data['Probability'] = data['bfactor']

    georep = ghm.GeoHTML(title, html_name, cols=3)
    georep.addLineComment(name + ', histograms')
    for geox in geos:
        print(geox)
        georep.addPlot1d(data,'histogram',geox)

    #georep.changeColNumber(len(pairs))
    print(' ...  Line 1 Pairs  #############################')
    georep.addLineComment(name + ', scatters')
    for geox,geoy in pairs:
        print(geox,geoy)
        georep.addPlot2d(data,'seaborn',geo_x=geox,geo_y=geoy,hue='aa')

    #filter extremes
    georep.changeColNumber(1)
    for geo in geos:
        print(' ...',geo,'extremes...')
        data = data.sort_values(by=geo, ascending=False)
        dataMax = data[:10]
        dataMin = data[-10:]
        georep.addLineComment('Extreme ' +  geo + ' values')
        dataMax = dataMax[['pdb_code',geo,'chain','rid','aa','Probability']]
        dataMin = dataMin[['pdb_code', geo, 'chain', 'rid', 'aa', 'Probability']]
        georep.addDataFrame(dataMin)
        georep.addDataFrame(dataMax)


    georep.printReport()
    print("Saved report to",html_name)

def createHtmlFromTriples(name,ext,data,triples):
    html_name = 'Html/' + name + "_" + ext + '_Scatters.html'
    title = name + ": Explore Geometry from " + name
    print(name, ' Creating HTML files #############################')
    data['Probability'] = data['bfactor']

    georep = ghm.GeoHTML(title, html_name, cols=3)
    print(' ...  Line 1 triples  #############################')
    georep.addLineComment(name + ' triples')
    for geoA,geoB,geoC in triples:
        print(geoA,geoB,geoC)
        georep.addPlot2d(data,'scatter',geo_x=geoA,geo_y=geoB,hue=geoC,palette='rainbow')
        georep.addPlot2d(data, 'scatter', geo_x=geoB, geo_y=geoC, hue=geoA, palette='rainbow')
        georep.addPlot2d(data, 'scatter', geo_x=geoC, geo_y=geoA, hue=geoB, palette='rainbow')

    georep.printReport()
    print("Saved report to",html_name)