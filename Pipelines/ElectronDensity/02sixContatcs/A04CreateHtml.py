'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import sys
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Helpers as help
import numpy as np
from LeucipPy import GeoHTML as ghm

#***************************************************************************************************************
def runCreateReport(df_leuci, df_geometry,HtmlStart,HtmlEnd, OutID,GeoAx,Geos,Overlay):
    from LeucipPy import LeuciFile as lef

    print('LeucipPipeline:02SixContacts - 3) Creating html report ----------------------------------------------------')
    idxs = df_leuci.index
    start = 0
    end = len(idxs)
    if HtmlStart !=-1:
        start = HtmlStart
    if HtmlEnd !=-1:
        end = HtmlEnd

    title = OutID
    fileNameD = 'Html/Den_' + OutID
    fileNameR = 'Html/Rad_' + OutID
    fileNameL = 'Html/Lap_' + OutID
    if HtmlStart !=-1:
        title += " From:" + str(HtmlStart+1)
        fileNameD += '_' + str(HtmlStart+1)
        fileNameR += '_' + str(HtmlStart + 1)
        fileNameL += '_' + str(HtmlStart + 1)
    if HtmlEnd !=-1:
        title += " To:" + str(HtmlEnd)
        fileNameD += '__' + str(HtmlEnd)
        fileNameR += '__' + str(HtmlEnd)
        fileNameL += '__' + str(HtmlEnd)
    fileNameD += '.html'
    fileNameR += '.html'
    fileNameL += '.html'

    georepD = ghm.GeoHTML(title,fileNameD, cols=5)
    georepR = ghm.GeoHTML(title, fileNameR, cols=5)
    georepL = ghm.GeoHTML(title, fileNameL, cols=5)

    # first add a general geometry report
    georepD.addLineComment('Six Contact Analysis in 5 Planes. View=Density')
    georepR.addLineComment('Six Contact Analysis in 5 Planes. View=Radiant')
    georepL.addLineComment('Six Contact Analysis in 5 Planes. View=Laplacian')

    letters = ['B','C','D','E','F']


    for count in range(start,end):
        letter = letters[count%5]
        Geox = Geos[count % 5]
        idx = idxs[count]
        set = df_leuci['SET'][idx].lower()
        id = df_leuci['ID'][idx].lower()
        pdb = df_leuci['pdb_code'][idx].lower()
        count += 1
        datapath = 'Output/' + set + '_' + id + '_' + pdb + '_SLICESFILE.csv'

        # get the information we need from 03csv
        idcode = 'ID_A' + letter
        cls = help.getInfo(df_geometry, 'CLASS', pdb, idcode,id.upper())
        res = help.getInfo(df_geometry, 'resolution', pdb, idcode, id.upper())
        atomC = help.getInfo(df_geometry, 'atomA', pdb, idcode,id.upper())
        atomL = help.getInfo(df_geometry, 'atomB', pdb, idcode,id.upper())
        atomP = help.getInfo(df_geometry, 'atomC', pdb, idcode,id.upper())
        disA = help.getInfo(df_geometry, GeoAx, pdb, idcode,id.upper())
        disB = help.getInfo(df_geometry, Geox, pdb, idcode,id.upper())
        Coords1 = help.getInfo(df_geometry, 'Coords0', pdb, idcode,id.upper())
        Coords2 = help.getInfo(df_geometry, 'CoordsA', pdb, idcode,id.upper())
        Coords3 = help.getInfo(df_geometry, 'Coords' + letter, pdb, idcode,id.upper())
        cx,cy,cz = float(Coords1.split(':')[0]),float(Coords1.split(':')[1]),float(Coords1.split(':')[2])
        lx,ly,lz = float(Coords2.split(':')[0]),float(Coords2.split(':')[1]),float(Coords2.split(':')[2])
        px,py,pz = float(Coords3.split(':')[0]),float(Coords3.split(':')[1]),float(Coords3.split(':')[2])
        from LeucipPy import GeoCalculate as calc
        #print(cx,cy,cz)
        #print(lx, ly, lz)
        #print(px, py, pz)
        try:
            angABC = calc.getAngle(lx,ly,lz,cx,cy,cz,px,py,pz)
        except:
            angABC = 0



        if atomC != '':
            print("LeuciPipelines:01Images: processing", datapath, count, '/', len(idxs))
            try:
                lfil = lef.LeuciFile(datapath)
                df_den = lfil.getDataFrame('DENSITYSLICE_0')
                mtx_den = lfil.getDataFrameAsMatrix(df_den,'Density')[0]
                df_rad = lfil.getDataFrame('RADIANTSLICE_0')
                mtx_rad = lfil.getDataFrameAsMatrix(df_rad, 'Radiant')[0]
                df_lap = lfil.getDataFrame('LAPLACIANSLICE_0')
                mtx_lap = lfil.getDataFrameAsMatrix(df_lap, 'Laplacian')[0]
                df_pos = lfil.getDataFrame('POSITIONSLICE_0')
                mtx_pos = lfil.placeDataFrameOnMatrix(df_pos, 'Position',np.zeros(mtx_rad.shape))

                box = '\n' + cls
                box += '\n' + 'Res=' + str(res)
                box += '\n' + id.upper()
                box += '\n' + GeoAx + '=' + str(round(disA,4))
                box += '\n' + Geox + '=' + str(round(disB,4))
                box += '\n' + 'Angle' + '=' + str(round(angABC, 2))



                #georep.addBoxComment(box)

                mtx_den = np.ma.masked_where(mtx_den < -100, mtx_den)
                mtx_rad = np.ma.masked_where(mtx_rad < -100, mtx_rad)
                mtx_lap = np.ma.masked_where(mtx_lap < -100, mtx_lap)
                mtx_pos = np.ma.masked_where(mtx_pos == 0,mtx_pos)

                #georep.addContours(mtx_den, pdb, colourbar=False, palette='magma_r', cap=1,overlay=True)

                georepD.addSurface(mtx_den, pdb + " Density" + box, palette='magma_r', overlay=True,cap=0.2)
                georepD.addSurface(mtx_pos, pdb + " Density" + box, alpha=1, palette="inferno_r", overlay=True)
                georepD.addContours(mtx_den, pdb + " Density" + box,style='lines', colourbar=False, palette='magma', levels=8,alpha=0.8,linewidth=0.5,cap=0.2)

                georepR.addSurface(mtx_rad, pdb + " Radiant" + box, palette='bone', overlay=True)
                georepR.addSurface(mtx_pos, pdb + " Radiant" + box, alpha=1, palette="inferno_r")

                georepL.addSurface(mtx_lap, pdb + " Laplacian" + box, palette='magma', overlay=True,cap=-0.5)
                georepL.addSurface(mtx_pos, pdb + " Laplacian" + box, alpha=1, palette="inferno_r", overlay=True)
                georepL.addContours(mtx_lap, pdb + " Laplacian" + box, style='lines', colourbar=False, palette='magma_r',levels=8, alpha=0.8, linewidth=0.5,cap=-0.5)
            except:
                print("!!! LeuciPipelines:error", datapath, count, '/', len(idxs))

    georepD.printReport()
    georepR.printReport()
    georepL.printReport()
    print('HTML printed to',fileNameD,fileNameR,fileNameL)
