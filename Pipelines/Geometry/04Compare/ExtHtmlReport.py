# !/usr/bin/python3
import sys

import pandas as pd

from Class01Html import PlotThread


def runOneHtml(ordered_geos_csv, correlation_csv, html_filename,geoA,start_count,end_count,aa_inc,aa_exc):
    #print('success')
    df_ordered = pd.read_csv(ordered_geos_csv)
    allgeos = df_ordered['geo'].values
    stats = df_ordered['stat'].values
    dic_ordered = {}
    for i in range(0, len(allgeos)):
        dic_ordered[allgeos[i]] = stats[i]

    #with open(html_filename, 'w') as file:
    #    file.write('whatever')

    df_geometry = pd.read_csv(correlation_csv)
    if aa_inc != '':
        df_geometry = df_geometry.query("aa == '" + aa_inc + "'")
    if aa_exc != '':
        df_geometry = df_geometry.query("aa != '" + aa_exc + "'")

    classRun = PlotThread(html_filename,'x',geoA,'bfactor',aa_inc,'',25,int(start_count),int(end_count),allgeos,df_geometry)
    classRun.run()




if __name__ == '__main__':
    globals()['runOneHtml'](sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])