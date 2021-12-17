'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import Helpers as help
import numpy as np
from LeucipPy import GeoHTML as ghm

#***************************************************************************************************************
def runCreateReport(df_geometry,OutID,GeoAx,GeoBx,GeoCx,HB):
    from LeucipPy import LeuciFile as lef

    print('LeucipPipeline:Geometry:01Geometry6- 2) Creating html report ----------------------------------------------------')

    title = OutID
    fileName = 'Html/' + OutID + '.html'
    georep = ghm.GeoHTML(title,fileName, cols=3)

    georep.addLineComment('Rotating hue')
    georep.addPlot2d(df_geometry, 'scatter', geo_x=GeoAx, geo_y=GeoBx, hue=GeoCx, palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'scatter', geo_x=GeoBx, geo_y=GeoCx, hue=GeoAx, palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'scatter', geo_x=GeoCx, geo_y=GeoAx, hue=GeoBx, palette='gist_ncar')

    georep.addLineComment('RID Gap as hue (of x-axis)')
    df_geometry = df_geometry.sort_values(by='ridGapA', ascending=False)
    georep.addPlot2d(df_geometry, 'scatter', geo_x=GeoAx, geo_y=GeoBx, hue='ridGapA', palette='Set1',crange=[1,6])
    df_geometry = df_geometry.sort_values(by='ridGapB', ascending=False)
    georep.addPlot2d(df_geometry, 'scatter', geo_x=GeoBx, geo_y=GeoCx, hue='ridGapB', palette='Set1',crange=[1,6])
    df_geometry = df_geometry.sort_values(by='ridGapC', ascending=False)
    georep.addPlot2d(df_geometry, 'scatter', geo_x=GeoCx, geo_y=GeoAx, hue='ridGapC', palette='Set1',crange=[1,6])

    georep.addLineComment('Aa as hue (of x-axis)')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='aaA', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoBx, geo_y=GeoCx, hue='aaB', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoCx, geo_y=GeoAx, hue='aaC', palette='gist_ncar')

    georep.addLineComment('Atom as hue (of x-axis)')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='atomA', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoBx, geo_y=GeoCx, hue='atomB', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoCx, geo_y=GeoAx, hue='atomC', palette='gist_ncar')

    georep.printReport()
    print('...HTML printed to',fileName)


def runCreateReportAnalysis(df_geometry, OutID, GeoAx, GeoBx, GeoCx, HB):
    print('LeucipPipeline:Geometry:01Geometry6- 2) Creating html analysis ----------------------------------------------------')

    title = OutID
    fileName = 'Html/' + OutID + '.html'
    georep = ghm.GeoHTML(title, fileName, cols=3)

    georep.addPlot1d(df_geometry,'histogram',GeoAx,hue='aaA',bins=50)
    georep.addPlot1d(df_geometry, 'histogram', GeoBx, hue='aaB',bins=50)
    georep.addPlot1d(df_geometry, 'histogram', GeoCx, hue='aaC',bins=50)

    georep.addSeries(df_geometry[GeoAx].describe(),GeoAx,True)
    georep.addSeries(df_geometry[GeoBx].describe(), GeoBx, True)
    georep.addSeries(df_geometry[GeoCx].describe(), GeoCx, True)

    georep.printReport()
    print('...HTML printed to', fileName)

#**********************************************************************************************************
def runCreateReportAnalysisMulti(df_geometry, OutID, Geos):
    print('LeucipPipeline:Geometry:01Geometry6- 2) Creating html analysis ----------------------------------------------------')

    title = OutID
    fileName = 'Html/' + OutID + '.html'
    georep = ghm.GeoHTML(title, fileName, cols=3)

    for i in range(0,len(Geos),3):
        if i < len(Geos):
            #print(df_geometry[Geos[i]])
            xrange = (min(df_geometry[Geos[i]].values),max(df_geometry[Geos[i]].values))
            georep.addPlot1d(df_geometry,'histogram',Geos[i],hue='pdb_code',bins=50,xrange=xrange)
        else:
            georep.addBoxComment('')
        if i+1 < len(Geos):
            xrange = (min(df_geometry[Geos[i+1]].values), max(df_geometry[Geos[i+1]].values))
            georep.addPlot1d(df_geometry,'histogram',Geos[i+1],hue='pdb_code',bins=50,xrange=xrange)
        else:
            georep.addBoxComment('')
        if i+2 < len(Geos):
            xrange = (min(df_geometry[Geos[i+2]].values), max(df_geometry[Geos[i+2]].values))
            georep.addPlot1d(df_geometry,'histogram',Geos[i+2],hue='pdb_code',bins=50,xrange=xrange)
        else:
            georep.addBoxComment('')

        if i < len(Geos):
            georep.addSeries(df_geometry[Geos[i]].describe(),Geos[i],True)
        else:
            georep.addBoxComment('')
        if i+1 < len(Geos):
            georep.addSeries(df_geometry[Geos[i+1]].describe(), Geos[i+1], True)
        else:
            georep.addBoxComment('')
        if i+2 < len(Geos):
            georep.addSeries(df_geometry[Geos[i+2]].describe(), Geos[i+2], True)
        else:
            georep.addBoxComment('')

    georep.printReport()
    print('...HTML printed to', fileName)
