'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
from LeucipPy import GeoHTML as ghm

#***************************************************************************************************************
def runCreateReport(df_geometry,OutID,GeoAx,GeoBx,GeoCx,GeoDx,GeoEx,GeoFx,HB):
    from LeucipPy import LeuciFile as lef

    print('LeucipPipeline:Geometry:01Geometry6- 2) Creating html report ----------------------------------------------------')

    title = OutID
    fileName = 'Html/' + OutID + '.html'
    georep = ghm.GeoHTML(title,fileName, cols=3)

    georep.addLineComment('Top 2 CC with original atom as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='aa', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='atom', palette='gist_ncar')
    df_geometryp = df_geometry.sort_values(by='ridGapA', ascending=False)
    georep.addPlot2d(df_geometryp, 'scatter', geo_x=GeoAx, geo_y=GeoBx, hue='ridGapA', palette='Set1',crange=[1,6])

    georep.addLineComment('Top 2 CC with first contact as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='aaA', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='atomA', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='CLASS', palette='gist_ncar')

    georep.addLineComment('Top 2 CC with second contact as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='aaB', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='atomB', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='CLASS', palette='gist_ncar')

    georep.addLineComment('Third contact as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoBx, geo_y=GeoCx, hue='aaC', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoBx, geo_y=GeoCx, hue='atomC', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoBx, geo_y=GeoCx, hue='CLASS', palette='gist_ncar')


    georep.addLineComment('Fourth contact as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoCx, geo_y=GeoDx, hue='aaD', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoCx, geo_y=GeoDx, hue='atomD', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoCx, geo_y=GeoDx, hue='CLASS', palette='gist_ncar')


    georep.addLineComment('Fifth contact as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoDx, geo_y=GeoEx, hue='aaE', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoDx, geo_y=GeoEx, hue='atomE', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoDx, geo_y=GeoEx, hue='CLASS', palette='gist_ncar')

    georep.addLineComment('Sixth contact as hue')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoEx, geo_y=GeoFx, hue='aaF', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoEx, geo_y=GeoFx, hue='atomF', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoEx, geo_y=GeoFx, hue='CLASS', palette='gist_ncar')

    georep.addLineComment('IDs of close contacts')
    georep.changeColNumber(2)

    georep.addBoxComment('The 3rd contacts <=' + str(HB[2]))
    georep.addBoxComment('The 4th contacts <=' + str(HB[3]))

    df3 = df_geometry.query('`' + GeoCx + '` <= ' + str(HB[2]))
    if len(df3.index) > 0:
        df3 = df3[['pdb_code','CLASS', 'InfoC',GeoCx]]
        georep.addDataFrame(df3)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[2]))

    df4 = df_geometry.query('`' + GeoDx + '` <= ' + str(HB[3]))
    if len(df4.index) > 0:
        df4a = df4[['pdb_code','CLASS', 'InfoD',GeoDx]]
        georep.addDataFrame(df4a)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[3]))

    georep.addBoxComment('The 5th contacts <=' + str(HB[4]))
    georep.addBoxComment('The 6th contacts <=' + str(HB[5]))

    df5 = df_geometry.query('`' + GeoEx + '` <= ' + str(HB[4]))
    if len(df5.index) > 0:
        df5a = df5[['pdb_code','CLASS', 'InfoE',GeoEx]]
        georep.addDataFrame(df5a)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[4]))

    df6 = df_geometry.query('`' + GeoFx + '` <= ' + str(HB[5]))
    if len(df6.index) > 0:
        df6a = df6[['pdb_code','CLASS', 'InfoF',GeoFx]]
        georep.addDataFrame(df6a)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[5]))

    georep.changeColNumber(1)
    georep.addBoxComment('More info on the 4th contacts <=' + str(HB[3]))
    if len(df4.index) > 0:
        df4b = df4[['pdb_code','CLASS', 'InfoA','InfoB','InfoC','InfoD']]
        georep.addDataFrame(df4b)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[3]))

    georep.addBoxComment('More info on the 5th contacts <=' + str(HB[4]))
    if len(df5.index) > 0:
        df5b = df5[['pdb_code', 'CLASS', 'InfoA', 'InfoB', 'InfoC', 'InfoD', 'InfoE']]
        georep.addDataFrame(df5b)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[4]))

    georep.addBoxComment('More info on the 6th contacts <=' + str(HB[5]))
    if len(df6.index) > 0:
        df6b = df6[['pdb_code', 'CLASS', 'InfoA', 'InfoB', 'InfoC', 'InfoD', 'InfoE', 'InfoF']]
        georep.addDataFrame(df6b)
    else:
        georep.addBoxComment('Nothing <=' + str(HB[5]))


    georep.printReport()
    print('...HTML printed to',fileName)


def runCreateReportAnalysis(df_geometry, OutID, GeoAx, GeoBx, GeoCx, GeoDx, GeoEx, GeoFx, HB):
    print('LeucipPipeline:Geometry:01Geometry6- 2) Creating html analysis ----------------------------------------------------')

    title = OutID
    fileName = 'Html/' + OutID + '.html'
    georep = ghm.GeoHTML(title, fileName, cols=3)

    #split them out
    df1 = df_geometry.query('`' + GeoAx + '` <= ' + str(HB[0]))
    df2 = df_geometry.query('`' + GeoBx + '` <= ' + str(HB[1]))
    df3 = df_geometry.query('`' + GeoCx + '` <= ' + str(HB[2]))
    df4 = df_geometry.query('`' + GeoDx + '` <= ' + str(HB[3]))
    df5 = df_geometry.query('`' + GeoEx + '` <= ' + str(HB[4]))
    df6 = df_geometry.query('`' + GeoFx + '` <= ' + str(HB[5]))


    total_num = len(df_geometry.index)
    if total_num > 0:
        num1 = len(df1.index)
        num2 = len(df2.index)
        num3 = len(df3.index)
        num4 = len(df4.index)
        num5 = len(df5.index)
        num6 = len(df6.index)
        num0 = total_num - num1

        withNoContacts = int(num0*100/total_num)
        with1Contact = int(100*(num1-num2) / total_num)
        with2Contacts = int(100*(num2 - num3) / total_num)
        with3Contacts = int(100*(num3 - num4) / total_num)
        with4Contacts = int(100*(num4 - num5) / total_num)
        with5Contacts = int(100*(num5 - num6) / total_num)
        with6Contacts = int(100*(num6) / total_num)

        georep.addPlot1d(df_geometry,'histogram',GeoAx,hue='aaA',bins=50)
        georep.addPlot1d(df_geometry, 'histogram', GeoBx, hue='aaB',bins=50)
        georep.addPlot1d(df_geometry, 'histogram', GeoCx, hue='aaC',bins=50)

        georep.addSeries(df_geometry[GeoAx].describe(),GeoAx,True)
        georep.addSeries(df_geometry[GeoBx].describe(), GeoBx, True)
        georep.addSeries(df_geometry[GeoCx].describe(), GeoCx, True)

        box = 'Having no Contacts <=' + str(HB[0]) + '<br/>' + str(num0) + ' / ' + str(total_num) + ' = ' + str(withNoContacts) + "%" + '<br/><br/>'
        box += 'Having 1 Contact <=' + str(HB[0]) + '<br/>' + str(num1-num2) + ' / ' + str(total_num) + ' = ' + str(with1Contact) + "%"
        georep.addBoxComment(box)
        georep.addBoxComment('Having 2 Contacts <=' + str(HB[1]) + '<br/>' + str(num2-num3) + ' / ' + str(total_num) + ' = ' + str(with2Contacts) + "%")
        georep.addBoxComment('Having 3 Contacts <=' + str(HB[2]) + '<br/>' + str(num3-num4) + ' / ' + str(total_num) + ' = ' + str(with3Contacts) + "%")


        georep.addPlot1d(df_geometry, 'histogram', GeoDx, hue='aaD',bins=50)
        georep.addPlot1d(df_geometry, 'histogram', GeoEx, hue='aaE',bins=50)
        georep.addPlot1d(df_geometry, 'histogram', GeoFx, hue='aaF',bins=50)

        georep.addSeries(df_geometry[GeoDx].describe(), GeoDx, True)
        georep.addSeries(df_geometry[GeoEx].describe(), GeoEx, True)
        georep.addSeries(df_geometry[GeoFx].describe(), GeoFx, True)

        georep.addBoxComment('Having 4 Contacts <=' + str(HB[3]) + '<br/>' + str(num4-num5) + ' / ' + str(total_num) + ' = ' + str(with4Contacts) + "%")
        georep.addBoxComment('Having 5 Contacts <=' + str(HB[4]) + '<br/>' + str(num5-num6) + ' / ' + str(total_num) + ' = ' + str(with5Contacts) + "%")
        georep.addBoxComment('Having 6 Contacts <=' + str(HB[5]) + '<br/>' + str(num6) + ' / ' + str(total_num) + ' = ' + str(with6Contacts) + "%")



    georep.printReport()
    print('...HTML printed to', fileName)
