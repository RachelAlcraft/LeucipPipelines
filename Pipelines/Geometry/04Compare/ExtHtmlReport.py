# !/usr/bin/python3
import sys

import pandas as pd

from Class01Html import PlotThread


def runOneHtml(ordered_geos_csv, correlation_csv, html_filename,geoA,start_count,end_count):
    #print('success')
    df_ordered = pd.read_csv(ordered_geos_csv)
    allgeos = df_ordered['geo'].values
    stats = df_ordered['stat'].values
    dic_ordered = {}
    for i in range(0, len(allgeos)):
        dic_ordered[allgeos[i]] = stats[i]

    with open(html_filename, 'w') as file:
        file.write('whatever')

    df_geometry = pd.read_csv(correlation_csv)

    classRun = PlotThread(html_filename,'x',geoA,'aa','ALA','',25,int(start_count),int(end_count),allgeos,df_geometry)
    classRun.run()




if __name__ == '__main__':
    globals()['runOneHtml'](sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])