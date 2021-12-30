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
    georep.changeColNumber(len(pairs))
    print(' ...  Line 1 Pairs  #############################')
    georep.addLineComment(name + ' Proteome, histograms')
    for geox,geoy in pairs:
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

