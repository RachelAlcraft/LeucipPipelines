'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import pandas as pd
from LeucipPy import GeoHTML as ghm
from LeucipPy import GeoDataFrame as gdf


#***************************************************************************************************************
def runCreateReport(ID, ID_high, ID_redo, df_geometry_high, df_geometry_redo, geosHist, triples):
    html_name = 'Html/' + ID + '_Compare.html'
    title = ID + ": Compare Geometry"
    print(ID, ' Creating HTML files #############################')

    georep = ghm.GeoHTML(title, html_name, cols=2)

    high_pdb = df_geometry_high['pdb_code'].unique()
    redo_pdb = df_geometry_redo['pdb_code'].unique()
    georep.addBoxComment('High has ' + str(len(high_pdb)) + ' pdbs')
    georep.addBoxComment('Redo has ' + str(len(redo_pdb)) + ' pdbs')

    print(' ...  Line 1 triples  #############################')
    for geo in geosHist:
        print('hist',geo)
        georep.addPlot1d(df_geometry_high,'histogram',geo,title=ID_high,bins=50,hue='pdb_code',xrange=[min(df_geometry_high[geo]),max(df_geometry_high[geo])])
        georep.addPlot1d(df_geometry_redo, 'histogram', geo, title=ID_redo,bins=50,hue='pdb_code',xrange=[min(df_geometry_redo[geo]),max(df_geometry_redo[geo])])
        georep.addSeries(df_geometry_high[geo].describe(), geo + ' ' + ID_high,True)
        georep.addSeries(df_geometry_redo[geo].describe(), geo + ' ' + ID_redo,True)

        high_summary = df_geometry_high.groupby('aa')[geo].describe()
        high_summary = high_summary.round(decimals=4)
        high_summary['aa'] =high_summary.index
        georep.addDataFrame(high_summary)

        redo_summary = df_geometry_redo.groupby('aa')[geo].describe()
        redo_summary = redo_summary.round(decimals=4)
        redo_summary['aa'] = redo_summary.index
        georep.addDataFrame(redo_summary)



    #georep.addLineComment(ID + ' triples')
    for geoA, geoB, geoC in triples:
        print('triple',geoA, geoB, geoC)
        if geoC in ['aa','HB_atom','HB_aa']:
            print('Triple seaborn',geoA, geoB, geoC)
            georep.addPlot2d(df_geometry_high, 'seaborn', geo_x=geoA, geo_y=geoB, hue=geoC, title=ID_high,palette='rainbow')
            georep.addPlot2d(df_geometry_redo, 'seaborn', geo_x=geoA, geo_y=geoB, hue=geoC, title=ID_redo,palette='rainbow')
        else:
            print('Triple scatter', geoA, geoB, geoC)
            georep.addPlot2d(df_geometry_high, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoC, title=ID_high, palette='rainbow')
            georep.addPlot2d(df_geometry_redo, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoC, title=ID_redo, palette='rainbow')


    georep.printReport()
    print("Saved report to", html_name)

def runCreateBoundariesReport(set_id, df_list, geosHist, triples):
    html_name = 'Html/' + set_id + '_Boundaries.html'
    title = set_id + ": Compare Geometry"
    print(set_id, ' Creating HTML files #############################')

    georep = ghm.GeoHTML(title, html_name, cols=len(df_list))

    georep.addLineComment('This report contains 3 columns for original data, and randomly sampled data of both psi and tau.')
    georep.addLineComment('The first reports are histograms, followed by scatter plots')
    georep.addLineComment('Histogram Reports for =' + str(geosHist))
    scat_dic = {}
    scat_dic['x'] = []
    scat_dic['y'] = []
    scat_dic['hue'] = []
    for geoA,geoB,geoC,l,u in triples:
        scat_dic['x'].append(geoA)
        scat_dic['y'].append(geoB)
        scat_dic['hue'].append(geoC)
        sct_df = pd.DataFrame.from_dict(scat_dic)

    georep.changeColNumber(1)
    georep.addDataFrame(sct_df, title='Scatter Reports')
    georep.addLineComment('')
    georep.changeColNumber(len(df_list))
    georep.addBoxComment('Original Data')
    georep.addBoxComment('Randomly sampled even PSI bands every 30 degrees')
    georep.addBoxComment('Randomly sampled even TAU bands every 2 degrees from 110-116')
    print(' ...  Line 1 triples  #############################')
    for geo in geosHist:
        print('hist',geo)
        georep.addLineComment('Histograms for ' + geo)
        for idg,df in df_list:
            georep.addPlot1d(df,'histogram',geo,title=idg,bins=50,hue='pdb_code',xrange=[min(df[geo]),max(df[geo])])
        for idg,df in df_list:
            georep.addSeries(df[geo].describe(), geo + ' ' + idg,True)

    #georep.addLineComment(ID + ' triples')
    for geoA, geoB, geoC, lower,upper in triples:
        georep.addLineComment('Scatters for ' + geoA + '-' + geoB + '-' + geoC)
        print('triple',geoA, geoB, geoC)
        for idg, df in df_list:
            df_use = df
            df_use = df.query("aa != 'PRO'")
            if geoC == 'PROB':
                georep.addPlot2d(df_use, 'probability', geo_x=geoA, geo_y=geoB, hue=geoC, title=idg,palette='cubehelix_r')
            elif geoC in ['aa','HB_atom','HB_atom2','HB_aa'] or geoB in ['aa','HB_atom','HB_atom2','HB_aa']:
                df_use = df_use.sort_values(by=geoC, ascending=True)
                georep.addPlot2d(df_use, 'seaborn', geo_x=geoA, geo_y=geoB, hue=geoC, title=idg,palette='tab20',crange=[lower,upper])
            elif geoC[0] == '@':
                geo = geoC[1:]
                georep.changeColNumber(5)
                for am in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']:
                    df_use = df.query("aa == '" + am + "'")
                    georep.addPlot2d(df_use, 'scatter', geo_x=geoA, geo_y=geoB, hue=geo, title=am + ' ' + idg, palette='rainbow',crange=[lower, upper])
                georep.changeColNumber(len(df_list))
            else:
                georep.addPlot2d(df_use, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoC, title=idg, palette='rainbow',crange=[lower,upper])






    georep.printReport()
    print("Saved report to", html_name)