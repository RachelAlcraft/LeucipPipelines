from threading import Thread
import random
import A04DifferenceImage as A04
import pandas as pd


#class PlotThread(Thread):
class PlotThread():
    def __init__(self,html_name,ID,geoA,hue,aa,tag,outlier_cut,start_count,end_count,allgeos,df_geometryB):
        #Thread.__init__(self)
        self.html_name = html_name

        self.ID = ID
        self.geoA = geoA
        self.hue = hue
        self.aa = aa
        self.tag = tag
        self.outlier_cut = outlier_cut
        self.start_count = start_count
        self.end_count = end_count
        self.allgeos = allgeos
        self.df_geometryB = df_geometryB

    def run(self):  # the processing function
        print('---Inside thread for',self.ID,self.geoA,self.start_count,self.end_count,len(self.allgeos))
        html_name = self.html_name

        ID = self.ID
        geoA = self.geoA
        hue = self.hue
        aa = self.aa
        tag = self.tag
        outlier_cut = self.outlier_cut
        start_count = self.start_count
        end_count = self.end_count
        allgeos = self.allgeos
        df_geometryB = self.df_geometryB

        geos = allgeos[start_count:end_count]
        list_by_stats = []

        count = 0
        for geoB in geos:
            count += 1
            if geoA != geoB:
                # First make the appropriate distributions
                #print('---Make distribution difference for', geoB, count, '/', len(geos))
                df_col = df_geometryB.sort_values(by=geoB, ascending=True)
                df_col = df_col.iloc[outlier_cut:, :]
                df_col = df_col.sort_values(by=geoB, ascending=False)
                df_col = df_col.iloc[outlier_cut:, :]
                arA, arB, arDiff, minv, maxv, stat, arDiffSq = A04.createIdealisedDifferencePlot(geoA, geoB, df_col)

                # create a random scatter just for visual check
                cutA = list(df_col[geoA].values)
                random.shuffle(cutA)
                cutB = list(df_col[geoB].values)
                random.shuffle(cutB)
                dic_cut = {}
                dic_cut[geoA] = cutA
                dic_cut[geoB] = cutB
                df_cut = pd.DataFrame.from_dict(dic_cut)
                df_cut = df_cut.sort_values(by=geoA, ascending=False)
                # arA1, arB1, arDiff1, minv1, maxv1, stat1 = A04.createDifferencePlot(geoA, geoB, df_col, df_cut)
                # sorted_list_by_stats.append([geoB, arA, arB, arDiff, minv, maxv, stat, arA1, arB1, arDiff1, minv1, maxv1, stat1])
                list_by_stats.append([geoB, arA, arB, arDiff, minv, maxv, stat, arDiffSq])

        #html_name = 'Html/' + ID + aa + tag + '_Dependency_' + str(start_count) + '_' + str(end_count) + '.html'
        title = geoA + ": search for dependent variables, using correlation metric"
        print('---',ID, ' Creating HTML files #############################')

        from LeucipPy import GeoHTML as ghm
        georep = ghm.GeoHTML(title, html_name, cols=5)

        count = 0
        # for geoB,arA, arB, arDiff, minv, maxv, stat,arA1, arB1, arDiff1, minv1, maxv1, stat1 in sorted_list_by_stats:
        for geoB, arA, arB, arDiff, minv, maxv, stat, arDiffSq in list_by_stats:
            count += 1

            print('---Making html', geoA, geoB, count, '/', end_count, '/', len(allgeos))

            georep.addLineComment('Geos= (' + geoA + ' , ' + geoB + ') - Correlation statistic=' + str(round(stat, 5)))

            df_col = df_geometryB.sort_values(by=geoB, ascending=True)
            df_col = df_col.iloc[outlier_cut:, :]
            df_col = df_col.sort_values(by=geoB, ascending=False)
            df_col = df_col.iloc[outlier_cut:, :]
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
                georep.addPlot2d(df_col, plottype='seaborn', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow',
                                 title=geoB + ' Original')
                georep.addPlot2d(df_cut, plottype='seaborn', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow',
                                 title=geoB + ' Randomised')
            else:
                georep.addPlot2d(df_col, plottype='scatter', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow',
                                 title=geoB + ' Original')
                georep.addPlot2d(df_cut, plottype='scatter', geo_x=geoA, geo_y=geoB, hue=hue, palette='rainbow',
                                 title=geoB + ' Randomised')

            georep.addSurface(arA, 'Orig plot', cmin=0, cmax=maxv, palette='Blues', colourbar=True)
            # georep.addSurface(arDiff1, 'Rand diff plot=' + str(round(stat1, 5)), cmin=minv1, cmax=maxv1, palette='seismic',colourbar=True)
            # georep.addSurface(arB1, 'Rand plot', cmin=0, cmax=maxv1, palette='Blues', colourbar=True)
            georep.addSurface(arDiff, 'Diff plot=' + str(round(stat, 5)), cmin=minv, cmax=maxv, palette='RdBu',
                              colourbar=True)
            georep.addSurface(arB, 'Convolved plot', cmin=0, cmax=maxv, palette='Reds', colourbar=True)


        print('---Saved to', html_name)
        georep.printReport()

