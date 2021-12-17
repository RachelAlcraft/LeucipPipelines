'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
LeucipPipelines:ElectronDensity:01Images
This pipeline takes
- a list of pdb structures
- some geoemtric parameters
- 3 atoms for a plane
- some inclusions or exclusions on the data set based on atoms, amino acids and geoemtric paramaters
It outputs
- An html report that contains
--the electron density
--radiant image for each sample
 --a geometry report
'''

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# We can run stages in case of failure recovery
runCreateData, runCreateElectronDensity, runCreateReport = True,True,True
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import Config as cfg
import Helpers as help
import subprocess as sub
import shutil
import Bio.PDB as bio
import pandas as pd
import numpy as np
from LeucipPy import GeoHTML as ghm
from LeucipPy import GeoDataFrame as gdf
from LeucipPy import LeucipPy as leu
from LeucipPy import GeoCalculate as calc



#***************************************************************************************************************
if runCreateData:
    print('LeucipPipeline:01Images - 1) Creating data -----------------------------------------------------------')
    pdb_df = pd.read_csv(cfg.PdbFile)
    pdb_df = pdb_df.sort_values(by=['PDB'], ascending=True)
    pdbs = pdb_df['PDB'].values
    ### 2) Loop through and create bio python objects
    parser = bio.PDBParser()
    strucs = []
    count = 1
    if cfg.CAP!=-1:#while testing reduce data
        pdbs = pdbs[:cfg.CAP]
    for pdb in pdbs:
        pdb = pdb.lower()
        print("CL: creating structures", count, '/', len(pdbs))
        count += 1
        pdb_file, pdb_html_loc = leu.getPdbLink(pdb)
        struc = parser.get_structure(pdb, cfg.PdbDirectory + pdb_file)
        strucs.append(struc)

    geo = gdf.GeoDataFrame(strucs, log=0)
    geos = [cfg.GeoA, cfg.GeoB]
    data = geo.calculateGeometry(geos, log=1)
    data.to_csv("Csv/" + cfg.ID + "_01_Geometry.csv", index=False)

    data[cfg.GeoAx] = data[cfg.GeoA]
    data[cfg.GeoBx] = data[cfg.GeoB]

    # make an id out of the associated atoms
    data['atom1'] = data.apply(lambda row: help.applyATOM(row['info' + cfg.GeoA], 0), axis=1)
    data['rid2'] = data.apply(lambda row: help.applyRID(row['info' + cfg.GeoA]), axis=1)
    data['aa2'] = data.apply(lambda row: help.applyAA(row['info' + cfg.GeoA]), axis=1)
    data['chain2'] = data.apply(lambda row: help.applyCHAIN(row['info' + cfg.GeoA]), axis=1)
    data['atom2'] = data.apply(lambda row: help.applyATOM(row['info' + cfg.GeoA],1), axis=1)

    data['rid3'] = data.apply(lambda row: help.applyRID(row['info' + cfg.GeoB]), axis=1)
    data['aa3'] = data.apply(lambda row: help.applyAA(row['info' + cfg.GeoB]), axis=1)
    data['chain3'] = data.apply(lambda row: help.applyCHAIN(row['info' + cfg.GeoB]), axis=1)
    data['atom3'] = data.apply(lambda row: help.applyATOM(row['info' + cfg.GeoB],1), axis=1)

    data['CLASS'] = data.apply(lambda row: help.applyCLASS(pdb_df,row['pdb_code']), axis=1)
    data['ID'] = data.apply(lambda row: help.applyID(row['aa'],row['rid'],row['chain'],row['atom1'],row['aa2'],row['rid2'],row['chain2'],row['atom2'],row['aa3'],row['rid3'],row['chain3'],row['atom3']), axis=1)

    data.to_csv("Csv/" + cfg.ID + "_02_Geometry.csv", index=False)

    if cfg.GeoAMin != -1:
        data = data.query('`' + cfg.GeoA + '` >= ' + str(cfg.GeoAMin))
    if cfg.GeoAMax != -1:
        data = data.query('`' + cfg.GeoA + '` <= ' + str(cfg.GeoAMax))
    if cfg.GeoBMin != -1:
        data = data.query('`' + cfg.GeoB + '` >= ' + str(cfg.GeoBMin))
    if cfg.GeoBMax != -1:
        data = data.query('`' + cfg.GeoB + '` <= ' + str(cfg.GeoBMax))
    if cfg.Inclusions != {}:
        data = geo.filterDataFrame(data, inclusions=cfg.Inclusions)
    if cfg.Exclusions != {}:
        data = geo.filterDataFrame(data, exclusions=cfg.Exclusions)

    atomsdata = geo.calculateData(log=1)
    data['Coords1'] = data.apply(lambda row: help.applyCOORDS(geo, atomsdata, row['pdb_code'], row['info' + cfg.GeoA], 0),axis=1)
    data['Coords2'] = data.apply(lambda row: help.applyCOORDS(geo, atomsdata, row['pdb_code'], row['info' + cfg.GeoA], 1),axis=1)
    data['Coords3'] = data.apply(lambda row: help.applyCOORDS(geo, atomsdata, row['pdb_code'], row['info' + cfg.GeoB], 1),axis=1)
    data[cfg.GeoCx] = data.apply(lambda row: help.applyDISTANCE(row['Coords2'], row['Coords3']), axis=1)

    if cfg.GeoCMin != -1:
        data = data.query('`' + cfg.GeoCx + '` >= ' + str(cfg.GeoCMin))
    if cfg.GeoCMax != -1:
        data = data.query('`' + cfg.GeoCx + '` <= ' + str(cfg.GeoCMax))

    data.to_csv("Csv/" + cfg.ID + "_03_Geometry.csv", index=False)
    # Now create the coords commands

    results_tup = []
    results_cols = ['SET','ID','pdb_code', 'CMD']
    idxs = data.index
    for idx in idxs:
        pdb_code = data['pdb_code'][idx].lower()
        chain1, aa1, rid1, atom1, coords1 = data['chain'][idx], data['aa'][idx], data['rid'][idx], data['atom1'][idx],data['Coords1'][idx]
        chain2, aa2, rid2, atom2, coords2 = data['chain2'][idx], data['aa2'][idx], data['rid2'][idx], data['atom2'][idx],data['Coords2'][idx]
        chain3, aa3, rid3, atom3,coords3 = data['chain3'][idx], data['aa3'][idx], data['rid3'][idx], data['atom3'][idx],data['Coords3'][idx]
        id1 = str(aa1) + str(rid1) + str(chain1) + "." + str(atom1) + "_" + str(aa2) + str(rid2) + str(chain2) + "." + str(atom2) + "_" + str(aa3) + str(rid3) + str(chain3) + "." + str(atom3)

        cx, cy, cz = float(coords1.split(':')[0]),float(coords1.split(':')[1]),float(coords1.split(':')[2])
        lx,ly,lz = float(coords2.split(':')[0]),float(coords2.split(':')[1]),float(coords2.split(':')[2])
        px, py, pz = float(coords3.split(':')[0]),float(coords3.split(':')[1]),float(coords3.split(':')[2])
        cx,cy,cz = round(cx,5),round(cy,5),round(cz,5)
        lx, ly, lz = round(lx, 5), round(ly, 5), round(lz, 5)
        px, py, pz = round(px, 5), round(py, 5), round(pz, 5)
        #"SLICESFILE|3nir|5|1|0|20_0.1|9.71_-12.376_15.907:0.975_-15.224_7.167:2.52_-13.584_7.06"
        CMD = "SLICESFILE|" + pdb_code  + "|5|1|0|15_0.075|"
        CMD +=  str(cx) + "_" + str(cy) + "_" + str(cz)
        CMD +=  ":" + str(lx) + "_" + str(ly) + "_" + str(lz)
        CMD += ":" + str(px) + "_" + str(py) + "_" + str(pz)

        results_tup.append([cfg.ID,id1,pdb_code,CMD])

    ### 7) Save the dataframe
    data_for_leuci = pd.DataFrame(results_tup,columns=['SET','ID','pdb_code','CMD'])
    data_for_leuci.to_csv("Csv/" + cfg.ID + "_04_LeucipPlus.csv", index=False)
    print("LeucipPipelines:01Images - Saved to", "Csv/" + cfg.ID + "_04_LeucipPlus.csv")


#***************************************************************************************************************
if runCreateElectronDensity:
    print('LeucipPipeline:01Images - 2) Creating electron density -----------------------------------------------')
    df_leuci = pd.read_csv("Csv/" + cfg.ID + "_04_LeucipPlus.csv")

    idxs = df_leuci.index
    count = 0
    for idx in idxs:
        count += 1
        print("LeuciPipelines:01Images: running LeucipPlus.exe", count, '/', len(idxs))
        set = df_leuci['SET'][idx].lower()
        id = df_leuci['ID'][idx].lower()
        pdb = df_leuci['pdb_code'][idx].lower()

        cmd = df_leuci['CMD'][idx]
        pigP = sub.Popen([cfg.ExePath, cmd], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()
        outputpath = cfg.LogDirectory + pdb + '_SLICESFILE.csv'
        movepath = 'Output/'+ set + '_' + id + '_' + pdb + '_SLICESFILE.csv'
        shutil.move(outputpath, movepath)




#***************************************************************************************************************
if runCreateReport:
    from LeucipPy import GeoDataFrame as gdf
    from LeucipPy import LeuciFile as lef

    print('LeucipPipeline:01Images - 3) Creating html report ----------------------------------------------------')
    df_leuci = pd.read_csv("Csv/" + cfg.ID + "_04_LeucipPlus.csv")
    df_geometry = pd.read_csv("Csv/" + cfg.ID + "_03_Geometry.csv")
    idxs = df_leuci.index
    start = 0
    end = len(idxs)
    if cfg.HtmlStart !=-1:
        start = cfg.HtmlStart
    if cfg.HtmlEnd !=-1:
        end = cfg.HtmlEnd

    title = cfg.ID
    fileName = 'Html/' + cfg.ID
    if cfg.HtmlStart !=-1:
        title += " From:" + str(cfg.HtmlStart)
        fileName += '_' + str(cfg.HtmlStart)
    if cfg.HtmlEnd !=-1:
        title += " To:" + str(cfg.HtmlEnd)
        fileName += '__' + str(cfg.HtmlEnd)
    fileName += '.html'

    georep = ghm.GeoHTML(title,fileName, cols=4)

    # first add a general geometry report
    georep.addLineComment('Geometric analysis of the data')
    georep.addPlot2d(df_geometry,'scatter',geo_x=cfg.GeoAx,geo_y=cfg.GeoBx,hue=cfg.GeoCx,palette='Spectral')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=cfg.GeoAx, geo_y=cfg.GeoBx, hue='aa2', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=cfg.GeoAx, geo_y=cfg.GeoBx, hue='atom2', palette='gist_ncar')
    georep.addPlot2d(df_geometry, 'seaborn', geo_x=cfg.GeoAx, geo_y=cfg.GeoBx, hue='aa3', palette='gist_ncar')

    georep.changeColNumber(3)
    georep.addLineComment('Electron density of each observation')

    overlay_density = np.zeros([2,2])
    overlay_radiant = np.empty([2,2])
    first_overlay = True



    for count in range(start,end):
        idx = idxs[count]
        set = df_leuci['SET'][idx].lower()
        id = df_leuci['ID'][idx].lower()
        pdb = df_leuci['pdb_code'][idx].lower()
        count += 1
        datapath = 'Output/' + set + '_' + id + '_' + pdb + '_SLICESFILE.csv'
        print("LeuciPipelines:01Images: processing", datapath, count, '/', len(idxs))
        lfil = lef.LeuciFile(datapath)
        df_den = lfil.getDataFrame('DENSITYSLICE_0')
        mtx_den = lfil.getDataFrameAsMatrix(df_den,'Density')[0]
        df_rad = lfil.getDataFrame('RADIANTSLICE_0')
        mtx_rad = lfil.getDataFrameAsMatrix(df_rad, 'Radiant')[0]
        df_pos = lfil.getDataFrame('POSITIONSLICE_0')
        mtx_pos = lfil.placeDataFrameOnMatrix(df_pos, 'Position',np.zeros(mtx_rad.shape))

        if cfg.Overlay:
            if first_overlay:
                overlay_density = mtx_den
                overlay_radiant = mtx_rad
                first_overlay = False
            else:
                overlay_density = lfil.placeDataFrameOnMatrix(df_den, 'Density',overlay_density )
                overlay_radiant = lfil.placeDataFrameOnMatrix(df_rad, 'Radiant', overlay_radiant)


        # get the information we need from 03csv
        cls = help.getInfo(df_geometry, 'CLASS', pdb, id)
        atomC = help.getInfo(df_geometry,'atom1',pdb,id)
        atomL = help.getInfo(df_geometry, 'atom2',pdb,id)
        atomP = help.getInfo(df_geometry, 'atom3',pdb,id)
        disA = help.getInfo(df_geometry, cfg.GeoA,pdb,id)
        disB = help.getInfo(df_geometry, cfg.GeoB,pdb,id)
        disC = help.getInfo(df_geometry, cfg.GeoCx,pdb,id)

        box = "<b>" + pdb + "</b><br/><br/>" + cls.upper() + '<br/><br/>' + id.upper()
        box += '<br/><br/>' + 'Central=' + atomC + " Linear=" + atomL + " Planar=" + atomP
        box += '<br/><br/>' + cfg.GeoAx + '=' + str(round(disA,4))
        box += '<br/>' + cfg.GeoBx + '=' + str(round(disB,4))
        box += '<br/>' + cfg.GeoCx + '=' + str(round(disC,4))

        georep.addBoxComment(box)

        mtx_den = np.ma.masked_where(mtx_den < -100, mtx_den)
        mtx_rad = np.ma.masked_where(mtx_rad < -100, mtx_rad)
        mtx_pos = np.ma.masked_where(mtx_pos == 0,mtx_pos)

        georep.addContours(mtx_den, pdb, colourbar=False, palette='magma_r', cap=1,overlay=True)
        georep.addContours(mtx_den, pdb,style='lines', colourbar=False, palette='magma', cap=1,overlay=False)
        georep.addSurface(mtx_rad, id.upper(), palette='bone', overlay=True)
        georep.addSurface(mtx_pos, id.upper(), alpha=1, palette="inferno_r")

    if cfg.Overlay:
        georep.addLineComment('Overlay of all the above')
        georep.addBoxComment('This contains an overlay of all the single observations<br/>Note they are not currently normalised')
        georep.addContours(overlay_density, 'Overlay', colourbar=False, palette='magma_r', cap=1, overlay=True)
        georep.addContours(overlay_density, 'Overlay', style='lines', colourbar=False, palette='magma', cap=1, overlay=False)
        georep.addSurface(overlay_radiant, 'Overlay', palette='bone')


    georep.printReport()
    print('HTML printed to',fileName)
