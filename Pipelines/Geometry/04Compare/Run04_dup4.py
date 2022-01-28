#@@@ PIPELINE PARMAMETERS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import math
import random

order_geo_list,make_html = False,True
runs = [] #ID,csv,geoA,aa_inc,aa_exc,hue,tag,over_geos
runs.append(['Correlation_TAU','PW_High_02_Geometry.csv','N:CA:C',['ALL'],['PRO'],'aa','',[]])
runs.append(['Correlation_TAUm1','PW_High_02_Geometry.csv','C-1:N:CA',['ALL'],['PRO'],'aa','',[]])
outlier_cut = 25
chunk = 40
geos_to_abs = ['CA:C:O:N+1','CA-1:C-1:N:CA','CA:C:N+1:CA+1']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Some filtering of the data
aa_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
aa_NO_gly = ['ALA','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import pandas as pd
import sys
from datetime import datetime
sys.path.append('C:/Dev/Github/LeucipPipelines/Pipelines/1Library')
import Helpers as hlp
import A01MakeCsv as A01
import A02CreateHtml as A02
import A03FilterAndRandomise as A03
import A04DifferenceImage as A04
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
start = datetime.now()
for ID,csv,geoA,aa_inc,aa_exc,hue,tag,over_geos in runs:

    print('Loading dataframe')
    df_geometryAll = pd.read_csv("Csv/" + csv)
    if aa != '':
        df_geometryAll = df_geometryAll.query("aa == '" + aa + "'")
    for exaa in exclude_aa:
        df_geometryAll = df_geometryAll.query("aa != '" + exaa + "'")

    for gabs in geos_to_abs:
        df_geometryAll[gabs] =abs(df_geometryAll[gabs])

    #cut off the extreme outliers
    df_geometryAll = df_geometryAll.sort_values(by=geoA, ascending=True)
    df_geometryAll = df_geometryAll.iloc[outlier_cut:,:]
    df_geometryAll = df_geometryAll.sort_values(by=geoA, ascending=False)
    df_geometryAll = df_geometryAll.iloc[outlier_cut:,:]

    cols = []
    for col in df_geometryAll.columns:
        if ':' in col and 'info' not in col:
            cols.append(col)
            print(col)

    df_geometryA = df_geometryAll[cols]
    cols.append(hue)
    df_geometryB = df_geometryAll[cols]

    if order_geo_list:
        geo_list = df_geometryA.columns
        ordered_list_by_stats = []
        count = 0
        for geoB in geo_list:
            count += 1
            if geoA != geoB:
                # First make the appropriate distributions
                print('Make distribution difference for', geoB, count, '/', len(df_geometryA.columns))
                df_col = df_geometryA.sort_values(by=geoB, ascending=True)
                df_col = df_col.iloc[outlier_cut:, :]
                df_col = df_col.sort_values(by=geoB, ascending=False)
                df_col = df_col.iloc[outlier_cut:, :]
                arA, arB, arDiff, minv, maxv, stat,arDiffSq = A04.createIdealisedDifferencePlot(geoA, geoB, df_col)
                ordered_list_by_stats.append([geoB, stat])

        print('sort list to most dependent at the top')
        ordered_list_by_stats = sorted(ordered_list_by_stats, key=lambda tup: tup[1], reverse=True)

        ordered_dic = {}
        ordered_dic['geo'] = []
        ordered_dic['stat'] = []
        for geoB, stat in ordered_list_by_stats:
            ordered_dic['geo'].append(geoB)
            ordered_dic['stat'].append(stat)

        df_ordered = pd.DataFrame.from_dict(ordered_dic)
        df_ordered.to_csv("Csv/OrderedGeos_" + ID + aa + ".csv", index=False)

    if make_html:
        df_ordered = pd.read_csv("Csv/OrderedGeos_" + ID + aa + ".csv")
        geos = df_ordered['geo'].values[start_count-1:end_count]
        stats = df_ordered['stat'].values[start_count-1:end_count]
        dic_ordered = {}
        for i in range(0,len(geos)):
            dic_ordered[geos[i]] = stats[i]

        if override_geos != []:
            geos = override_geos
        list_by_stats = []

        count = 0
        for geoB in geos:
            count += 1
            if geoA != geoB:
                # First make the appropriate distributions
                print('Make distribution difference for', geoB, count, '/', len(geos))
                df_col = df_geometryB.sort_values(by=geoB, ascending=True)
                df_col = df_col.iloc[outlier_cut:, :]
                df_col = df_col.sort_values(by=geoB, ascending=False)
                df_col = df_col.iloc[outlier_cut:, :]
                arA, arB, arDiff, minv, maxv, stat, arDiffSq = A04.createIdealisedDifferencePlot(geoA, geoB, df_col)

                #create a random scatter just for visual check
                cutA = list(df_col[geoA].values)
                random.shuffle(cutA)
                cutB = list(df_col[geoB].values)
                random.shuffle(cutB)
                dic_cut = {}
                dic_cut[geoA] = cutA
                dic_cut[geoB] = cutB
                df_cut = pd.DataFrame.from_dict(dic_cut)
                df_cut = df_cut.sort_values(by=geoA, ascending=False)
                #arA1, arB1, arDiff1, minv1, maxv1, stat1 = A04.createDifferencePlot(geoA, geoB, df_col, df_cut)
                #sorted_list_by_stats.append([geoB, arA, arB, arDiff, minv, maxv, stat, arA1, arB1, arDiff1, minv1, maxv1, stat1])
                list_by_stats.append([geoB, arA, arB, arDiff, minv, maxv, stat,arDiffSq])

        print('sort list to most dependent at the top')

        #if override_geos == []:
        #    list_by_stats = sorted(list_by_stats, key=lambda tup: tup[6], reverse=True)


        html_name = 'Html/' + ID + aa + tag + '_Dependency_' + str(start_count) + '_' + str(end_count) + '.html'
        title = geoA + ": search for dependent variables, using correlation metric"
        print(ID, ' Creating HTML files #############################')

        from LeucipPy import GeoHTML as ghm
        georep = ghm.GeoHTML(title, html_name, cols=5)


        count = 0
        #for geoB,arA, arB, arDiff, minv, maxv, stat,arA1, arB1, arDiff1, minv1, maxv1, stat1 in sorted_list_by_stats:
        for geoB,arA, arB, arDiff, minv, maxv, stat,arDiffSq in list_by_stats:
            count +=1

            print('Making html',count, geoA, geoB, len(list_by_stats))
            georep.addLineComment('Geos= (' + geoA + ' , ' + geoB + ') - Correlation statistic=' + str(round(stat,5)))

            df_col = df_geometryB.sort_values(by=geoB, ascending=True)
            df_col = df_col.iloc[outlier_cut:,:]
            df_col = df_col.sort_values(by=geoB, ascending=False)
            df_col = df_col.iloc[outlier_cut:,:]
            df_col = df_col.sort_values(by=geoA, ascending=True)

            cutA = list(df_col[geoA].values)
            random.shuffle(cutA)
            cutB = list(df_col[geoB].values)
            random.shuffle(cutB)
            cutHue = list(df_col[hue].values)
            dic_cut = {}
            dic_cut[geoA] = cutA
            dic_cut[geoB] = cutB
            dic_cut[hue] = cutHue
            df_cut = pd.DataFrame.from_dict(dic_cut)
            df_cut = df_cut.sort_values(by=geoA, ascending=True)

            hues = df_cut.sort_values(by=hue, ascending=True)[hue].values

            if hue == 'aa':
                georep.addPlot2d(df_col, plottype='seaborn', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow',title=geoB+' Original')
                georep.addPlot2d(df_cut, plottype='seaborn', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow',title=geoB+' Randomised')
            else:
                georep.addPlot2d(df_col, plottype='scatter', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow', title=geoB + ' Original')
                georep.addPlot2d(df_cut, plottype='scatter', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow', title=geoB + ' Randomised')

            #arrrA = plt.hist2d(x, y, bins=10, range=None, density=False, weights=None,

            georep.addSurface(arA, 'Orig plot', cmin=0, cmax=maxv, palette='Blues',colourbar=True)
            #georep.addSurface(arDiff1, 'Rand diff plot=' + str(round(stat1, 5)), cmin=minv1, cmax=maxv1, palette='seismic',colourbar=True)
            #georep.addSurface(arB1, 'Rand plot', cmin=0, cmax=maxv1, palette='Blues', colourbar=True)
            georep.addSurface(arDiff, 'Diff plot=' + str(round(stat, 5)), cmin=minv, cmax=maxv, palette='RdBu', colourbar=True)
            georep.addSurface(arB, 'Convolved plot', cmin=0, cmax=maxv, palette='Reds', colourbar=True)

            #georep.addSurface(arDiffSq, 'Abs(Diff plot)=' + str(round(stat, 5)), cmin=0,cmax=1,palette='Greys',colourbar=True)


        print('Saved to',html_name)
        georep.printReport()

    end = datetime.now()
    hlp.printTime(start,end)


