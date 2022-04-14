'''
Rachel Alcraft: 09/02/2022
'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# INPUTS #
#control which steps you want to run
import numpy as np
modify_csv,recreate_html = True,True

runs = []
#runs.append(['Redo_GLY','hb'])
#runs.append(['Redo','hb'])
#runs.append(['High','hb'])
runs.append(['High_GLY','hb'])
#runs.append(['Redo_GLY_IDEAL','hb'])
runs.append(['High_GLY_IDEAL_PAIRED','hb'])
runs.append(['High_GLY_IDEAL_UNPAIRED','hb'])

for tag,tag2 in runs:
    csv_final = "Csv/PW_" + tag + "_01_Geometry.csv"
    html_filename = 'Html/11_HeatMap_' +tag+tag2+ '.html'
    title = 'Williams Divergence from Trivial: Heatmaps ' + tag + tag2 # used at the header of the html report
    csv_correlations = "Csv/10_DivCorr_" + tag + ".csv"
    csv_summary = "Csv/10_SumCorr_" + tag + ".csv"
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    from LeucipPy import HtmlReportMaker as hrm
    from LeucipPy import WilliamsDivergenceMaker as wcm
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import A0_Globals as glob
    ###############################################################################################
    data = pd.read_csv(csv_final)
    geos = []
    for col in data.columns:
        #if ':' in col and 'info' not in col and '(' not in col and '{' not in col:
        if ':' in col and 'info' not in col and '{' not in col:
        #if ':' in col and 'info' not in col:
            geos.append(col)
    if modify_csv:
        print('### Applying filters to csv data')
        data = data.query('occupancy == 1')
        data = data.query('bfactor <= 10')
        data = glob.trimDihs(data,15)

    if recreate_html:
        print('### Creating html reports')
        rep_mak = hrm.HtmlReportMaker(title,html_filename, cols=2)
        print('Create Williams Coefficient Maker')
        wcc = wcm.WilliamsDivergenceMaker(data,geos,density=1,log=1,norm=True,pval_iters=0)

        complete = wcc.getCoefficientsDataFrame()
        complete = complete.sort_values(by='stat', ascending=False)
        #complete.to_csv('Csv/10_Correlations.csv')
        summary = complete.groupby(by='geoA').sum()
        summary = summary.sort_values(by='stat')
        #summary.to_csv('Csv/10_Summary.csv')

        #complete = pd.read_csv(csv_correlations)
        #complete = complete.sort_values(by='stat', ascending=False)
        #summary = pd.read_csv(csv_summary)


        # remove the low correlations
        rem_df = complete[['geoA','geoB','stat']]
        ccut_off = 30
        len_cut_off = len(summary.index)-ccut_off
        for g in range(0,len_cut_off):
            geo = summary.index[g]
            rem_df = rem_df[rem_df['geoA'] != geo]
            rem_df = rem_df[rem_df['geoB'] != geo]

        red_piv = rem_df.pivot('geoA','geoB','stat')
        print(red_piv)

        for i in range(0,len(red_piv.columns)):
            red_piv = red_piv.sort_values(by=red_piv.columns[i], ascending=False)
            red_piv = red_piv[list(red_piv.index.values)]
        for i in range(len(red_piv.columns)-1,-1,-1):
            red_piv = red_piv.sort_values(by=red_piv.columns[i], ascending=False)
            red_piv = red_piv[list(red_piv.index.values)]

        fig, ax = plt.subplots()
        sns.heatmap(red_piv, annot=False, fmt='.2f', linewidth=.5, ax=ax, cmap='inferno_r', mask=red_piv.isnull(), vmin=0, xticklabels=red_piv.columns,yticklabels=red_piv.index)
        plt.title('Heatmap with least correlated removed')
        ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6, rotation=90)
        ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=6, rotation=0)
        ax.tick_params(which="both", bottom=True)
        rep_mak.addPlotOnly(fig, ax)

        #for frac,asc,ttl in [[1,False,'Heatmap of all geos'],[0.2,False,'Heatmap top 20%'],[0.2,True,'Heatmap bottom 20%']]:
        heat=wcc.getCoefficientsDataFrame(as_pivot=True,fraction=1,asc=False)

        for i in range(0, len(heat.columns)):
            heat = heat.sort_values(by=heat.columns[i], ascending=False)
            heat = heat[list(heat.index.values)]
        for i in range(len(heat.columns) - 1, -1, -1):
            heat = heat.sort_values(by=heat.columns[i], ascending=False)
            heat = heat[list(heat.index.values)]


        fig, ax = plt.subplots()
        sns.heatmap(heat, annot=False, fmt='.2f', linewidth=.5, ax=ax, cmap='inferno_r',mask=heat.isnull(),vmin=0)
        plt.title('Heatmap of all geos')
        ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6,rotation=90)
        ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=6, rotation=0)
        rep_mak.addPlotOnly(fig, ax)

        rep_mak.printReport()




