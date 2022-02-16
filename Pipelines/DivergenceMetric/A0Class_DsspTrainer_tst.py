import pandas as pd
import seaborn as sns

import A0Class_DsspTrainer as dssp
# check new data that we have not trained on
from LeucipPy import BioPythonMaker as bpm
from LeucipPy import DataFrameMaker as dfm
from LeucipPy import HtmlReportMaker as hrm

geos_for_report = ['N:O','CA-1:CA:CA+1','CA-1:CA+1','O:CA+1','CA:O-1','C:O','C:N+1','N:CA:C:N+1','C-1:N:CA:C']
aa = ''
recreate = True
train_geosAll=['N:O+2','N:O+3','N:O+4','O:N+2','O:N+3','O:N+4','O:{N}+2','N:{O}+2','N:CA:C','N:CA:C:N+1','C-1:N:CA:C']
train_geosA = ['N:CA:C:N+1','C-1:N:CA:C']
train_geosB=['O:{N}+2','N:{O}+2','N:CA:C','N:CA:C:N+1','C-1:N:CA:C']
train_geosC=['N:O+2','N:O+3','N:O+4','O:N+2','O:N+3','O:N+4']

outer_runs = []
#outer_runs.append([0.5,train_geosAll,'all'])
#outer_runs.append([0.5,train_geosA,'rama'])
#outer_runs.append([0.5,train_geosB,'hb2'])
outer_runs.append([0.5,train_geosC,'hb3'])
#outer_runs.append([0.3,train_geosB,'ca'])
#outer_runs.append([0.3,train_geosC,'four'])
#outer_runs.append([0.3,train_geosD,'five'])
#outer_runs.append([0.3,train_geosE,'all'])
#outer_runs.append([0.3,train_geosF,'co'])
#outer_runs.append([0.3,train_geosG,'omega'])

inner_runs = []
inner_runs.append(['QuickTest','C:/Dev/Github/ProteinDataFiles/pdb_traintest/','ent','pdb',0])
#inner_runs.append(['AlphaEcoli','C:/Dev/Github/ProteinDataFiles/pdb_alpha/EColi/','pdb','',50])
#inner_runs.append(['AlphaHuman','C:/Dev/Github/ProteinDataFiles/pdb_alpha/Human/','pdb','',50])

for fraction,training_geos,outer_tag in outer_runs:
    print(fraction, training_geos,outer_tag)
    report_geos = []
    for geo in geos_for_report:
        report_geos.append(geo)
    for geo in training_geos:
        if geo not in report_geos:
            report_geos.append(geo)

    # First train the classifier
    if recreate:
        trainingset = pd.read_csv('Csv/SSTrainingSet_hb.csv')
        print('### TRAINING THE CLASSIFIER ###')
        if training_geos == [] and aa == '':
            trainer = dssp.DsspTrainer(trainingset,fraction=fraction,log=1)
        elif aa == '':
            trainer = dssp.DsspTrainer(trainingset,fraction=fraction, log=1, geos=training_geos)
        else:
            trainer = dssp.DsspTrainer(trainingset,fraction=fraction, log=1, geos=training_geos,aa_list=[aa])

    for tag,dir,ext,pfx,count in inner_runs:
        print('#### running for',tag,dir,'##################')
        csv_final = "Csv/A0_TRAINING_Geometry_dssp_" + tag + str(fraction) + outer_tag + aa + ".csv"
        title = "Secondary Structure Classifier Test: " + tag+ str(fraction) + ' <br/>geos trained = '
        for geo in training_geos:
            title += str(geo) + ' , '
        html_filename = 'Html/A0_SecondaryTest_' + tag + str(fraction) +outer_tag +aa + '.html'

        if recreate:
            print('### Load structures from BioPython #############')
            strucs = bpm.loadPdbStructures([],dir,extension=ext,prefix=pfx,log=2,count=count)

            print('### Creating dataframe for correlations #############')
            geo_mak = dfm.DataFrameMaker(strucs, log=1)
            data = geo_mak.calculateGeometry(report_geos, log=1)

            if aa != '':
                data = data[data['aa'] == aa]

            print('##### ADD DSSP #####')
            data_dssp = trainer.classifySecondaryStructures(data)
            print(data_dssp)

            print('### Save dataframe ###############################')
            data_dssp.to_csv(csv_final, index=False)
            print("Saved to", csv_final)
        else:
            data_dssp = pd.read_csv(csv_final)

        print('### Creating html reports ########################')
        pal = sns.color_palette("bright")
        pal2 = []
        pal2.append('Gainsboro')
        for pl in pal:
            pal2.append(pl)
        customPalette = sns.set_palette(sns.color_palette(pal2))
        rep_mak = hrm.HtmlReportMaker(title, html_filename, cols=2)
        # rep_mak.addPlot1d(data_dssp,'histogram',geo_x='ss',hue='pdb_code',title='Secondary Structure Counts')
        rep_mak.addPlot2d(data_dssp, 'seaborn', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='ss', title='', palette=customPalette)
        rep_mak.addPlot2d(data_dssp, 'seaborn', 'N:CA:C:N+1', 'N:O', hue='ss', title='', palette=customPalette)
        rep_mak.addPlot2d(data_dssp, 'seaborn', 'C:O', 'C:N+1', hue='ss', title='', palette=customPalette)
        # rep_mak.addLineComment('Data by secondary structure')
        # rep_mak.changeColNumber(1)
        datagrouped = data_dssp[['ss', 'aa']]
        datagrouped = datagrouped.groupby(by=['ss']).count()
        datagrouped['ss'] = datagrouped.index
        datagrouped['count'] = datagrouped['aa']
        rep_mak.addDataFrame(datagrouped[['ss', 'count']])
        rep_mak.printReport()
        print('Exported html to', html_filename)

