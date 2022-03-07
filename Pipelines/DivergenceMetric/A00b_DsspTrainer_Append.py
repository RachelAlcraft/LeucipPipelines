import pandas as pd
import seaborn as sns
import A0Class_DsspTrainer as dssp

fraction = 0.5
recreate = True

if recreate:
    print('### training classifier ########################')
    trainingset = pd.read_csv('Csv/SSTrainingSet_hb.csv')
    trainer = dssp.DsspTrainer(trainingset,fraction=fraction)

#csvs = ["Redo_GLY","High_GLY","High","Redo",'SYN_GLY']
#csvs = ["Redo_GLY","High_GLY",'High_GLY_IDEAL','Redo_GLY_IDEAL']
csvs = ['High','High_GLY','High_IDEAL_PAIRED','High_IDEAL_UNPAIRED','High_GLY_IDEAL_PAIRED','High_GLY_IDEAL_UNPAIRED']

#csvs = ['SYN_GLY']
for tag in csvs:
    print('###', tag,'########################')
    csv = "Csv/PW_" + tag + "_01_Geometry.csv"
    html_filename = "Html/A0_TestClassifier_" + tag + str(fraction) + ".html"
    dssp_csv = csv + '_dssp.csv'
    print('reading', csv)
    data = pd.read_csv(csv)
    if recreate:
        print('### Creating dataframe ########################')
        data_dssp = trainer.classifySecondaryStructures(data)
        print(data_dssp)
        print('Save to', dssp_csv)
        data_dssp.to_csv(dssp_csv, index=False)
    else:
        data_dssp = pd.read_csv(dssp_csv)

    from LeucipPy import HtmlReportMaker as hrm
    print('### Creating html reports ########################')

    pal = sns.color_palette("bright")
    pal2 = []
    pal2.append('Gainsboro')
    for pl in pal:
        pal2.append(pl)
    customPalette = sns.set_palette(sns.color_palette(pal2))

    rep_mak = hrm.HtmlReportMaker("Test Classifier " + tag,html_filename, cols=4)
    rep_mak.addLineComment('New data by secondary structure')
    #rep_mak.addPlot1d(data_dssp,'histogram',geo_x='ss',hue='pdb_code',title='Secondary Structure Counts')
    rep_mak.addPlot2d(data_dssp,'seaborn','C-1:N:CA:C','N:CA:C:N+1',hue='ss',title='',palette=customPalette)
    rep_mak.addPlot2d(data_dssp,'seaborn','N:CA:C:N+1','N:O',hue='ss',title='',palette=customPalette)
    rep_mak.addPlot2d(data_dssp,'seaborn','C:O','C:N+1',hue='ss',title='',palette=customPalette)
    datagrouped = data_dssp[['ss','aa']]
    datagrouped = datagrouped.groupby(by=['ss']).count()
    datagrouped['ss'] = datagrouped.index
    datagrouped['count'] =datagrouped['aa']
    rep_mak.addDataFrame(datagrouped[['ss','count']])

    rep_mak.addLineComment('Prev data by aa')
    rep_mak.addPlot2d(data, 'seaborn', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='aa', title='', palette='Spectral_r')
    rep_mak.addPlot2d(data, 'seaborn', 'N:CA:C:N+1', 'N:O', hue='aa', title='', palette='Spectral_r')
    rep_mak.addPlot2d(data, 'seaborn', 'C:O', 'C:N+1', hue='aa', title='', palette='Spectral_r')
    rep_mak.addBoxComment('')

    rep_mak.printReport()
    print('Exported html to', html_filename)



