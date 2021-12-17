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
def runCreateReport(df_leuci, df_geometry,HtmlStart,HtmlEnd, OutID,GeoAx,GeoBx,GeoCx,Overlay):
    from LeucipPy import LeuciFile as lef

    print('LeucipPipeline:01Images - 3) Creating html report ----------------------------------------------------')
    idxs = df_leuci.index
    start = 0
    end = len(idxs)
    if HtmlStart !=-1:
        start = HtmlStart
    if HtmlEnd !=-1:
        end = HtmlEnd

    title = OutID
    fileName = 'Html/' + OutID
    if HtmlStart !=-1:
        title += " From:" + str(HtmlStart+1)
        fileName += '_' + str(HtmlStart+1)
    if HtmlEnd !=-1:
        title += " To:" + str(HtmlEnd)
        fileName += '__' + str(HtmlEnd)
    fileName += '.html'

    georep = ghm.GeoHTML(title,fileName, cols=4)

    # first add a general geometry report
    georep.addLineComment('Geometric analysis of the data')
    georep.addPlot2d(df_geometry,'scatter',geo_x=GeoAx,geo_y=GeoBx,hue=GeoCx,palette='Spectral')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='aa2', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='atom2', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=GeoAx, geo_y=GeoBx, hue='aa3', palette='gist_ncar')

    georep.changeColNumber(4)
    georep.addLineComment('Electron density of each observation')

    overlay_density = np.zeros([2,2])
    overlay_radiant = np.empty([2,2])
    overlay_laplacian = np.empty([2, 2])
    first_overlay = True
    overlay_msg = ''
    overlay_count = 0

    for count in range(start,end):
        idx = idxs[count]
        set = df_leuci['SET'][idx].lower()
        id = df_leuci['ID'][idx].lower()
        pdb = df_leuci['pdb_code'][idx].lower()
        count += 1
        datapath = 'Output/' + set + '_' + id + '_' + pdb + '_SLICESFILE.csv'

        # get the information we need from 03csv
        cls = help.getInfo(df_geometry, 'CLASS', pdb, id.upper())
        atomC = help.getInfo(df_geometry, 'atom1', pdb, id.upper())
        atomL = help.getInfo(df_geometry, 'atom2', pdb, id.upper())
        atomP = help.getInfo(df_geometry, 'atom3', pdb, id.upper())
        disA = help.getInfo(df_geometry, GeoAx, pdb, id.upper())
        disB = help.getInfo(df_geometry, GeoBx, pdb, id.upper())
        disC = help.getInfo(df_geometry, GeoCx, pdb, id.upper())
        Coords1 = help.getInfo(df_geometry, 'Coords1', pdb, id.upper())
        Coords2 = help.getInfo(df_geometry, 'Coords2', pdb, id.upper())
        Coords3 = help.getInfo(df_geometry, 'Coords3', pdb, id.upper())
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
            lfil = lef.LeuciFile(datapath)
            df_den = lfil.getDataFrame('DENSITYSLICE_0')
            mtx_den = lfil.getDataFrameAsMatrix(df_den,'Density')[0]
            df_rad = lfil.getDataFrame('RADIANTSLICE_0')
            mtx_rad = lfil.getDataFrameAsMatrix(df_rad, 'Radiant')[0]
            df_lap = lfil.getDataFrame('LAPLACIANSLICE_0')
            mtx_lap = lfil.getDataFrameAsMatrix(df_lap, 'Laplacian')[0]
            df_pos = lfil.getDataFrame('POSITIONSLICE_0')
            mtx_pos = lfil.placeDataFrameOnMatrix(df_pos, 'Position',np.zeros(mtx_rad.shape))

            if Overlay:
                if first_overlay:
                    overlay_density = mtx_den
                    overlay_radiant = mtx_rad
                    overlay_laplacian = mtx_lap
                    first_overlay = False
                else:
                    overlay_density = lfil.placeDataFrameOnMatrix(df_den, 'Density',overlay_density )
                    overlay_radiant = lfil.placeDataFrameOnMatrix(df_rad, 'Radiant', overlay_radiant)
                    overlay_laplacian = lfil.placeDataFrameOnMatrix(df_lap, 'Laplacian', overlay_laplacian)
                overlay_msg += '<br/>' + pdb + ':' + id
                overlay_count += 1




            box = "<b>" + pdb + "</b><br/><br/>" + cls.upper() + '<br/><br/>' + id.upper()
            box += '<br/><br/>' + 'Central=' + atomC + " Linear=" + atomL + " Planar=" + atomP
            box += '<br/><br/>' + GeoAx + '=' + str(round(disA,4))
            box += '<br/>' + GeoBx + '=' + str(round(disB,4))
            box += '<br/>' + GeoCx + '=' + str(round(disC,4))
            box += '<br/>' + 'Angle' + '=' + str(round(angABC, 2))
            box += '<br/><br/>' + 'Central' + '=' + str(Coords1)
            box += '<br/>' + 'Linear' + '=' + str(Coords2)
            box += '<br/>' + 'Planar' + '=' + str(Coords3)


            georep.addBoxComment(box)

            mtx_den = np.ma.masked_where(mtx_den < -100, mtx_den)
            mtx_rad = np.ma.masked_where(mtx_rad < -100, mtx_rad)
            mtx_lap = np.ma.masked_where(mtx_lap < -100, mtx_lap)
            mtx_pos = np.ma.masked_where(mtx_pos == 0,mtx_pos)

            #georep.addContours(mtx_den, pdb, colourbar=False, palette='magma_r', cap=1,overlay=True)
            georep.addSurface(mtx_den, pdb + " Density", palette='magma_r', overlay=True)
            georep.addSurface(mtx_pos, pdb + " Density", alpha=1, palette="inferno_r", overlay=True)
            georep.addContours(mtx_den, pdb + " Density",style='lines', colourbar=False, palette='magma', levels=8,alpha=0.8,linewidth=0.5)
            georep.addSurface(mtx_rad, pdb + " Radiant", palette='bone', overlay=True)
            georep.addSurface(mtx_pos, pdb + " Radiant", alpha=1, palette="inferno_r")
            georep.addSurface(mtx_lap, pdb + " Laplacian", palette='magma', overlay=True)
            georep.addSurface(mtx_pos, pdb + " Laplacian", alpha=1, palette="inferno_r", overlay=True)
            georep.addContours(mtx_lap, pdb + " Laplacian", style='lines', colourbar=False, palette='magma_r',levels=8, alpha=0.8, linewidth=0.5)

    if Overlay:
        box = 'This contains an overlay of all the single observations<br/>Note they are not currently normalised'
        box += '<br/> Obsrvations=' + str(overlay_count) + overlay_msg
        georep.addLineComment('Overlay of all the above')
        georep.addBoxComment(box)
        georep.addContours(overlay_density, 'Overlay Density', colourbar=False, palette='magma_r', cap=1, overlay=True)
        georep.addContours(overlay_density, 'Overlay Density', style='lines', colourbar=False, palette='magma', cap=1, overlay=False)
        georep.addSurface(overlay_radiant, 'Overlay Radiant', palette='bone')
        georep.addContours(overlay_laplacian, 'Overlay Laplacian', colourbar=False, palette='magma', cap=1, overlay=True)
        georep.addContours(overlay_laplacian, 'Overlay Laplacian', style='lines', colourbar=False, palette='magma_r')


    georep.printReport()
    print('HTML printed to',fileName)
