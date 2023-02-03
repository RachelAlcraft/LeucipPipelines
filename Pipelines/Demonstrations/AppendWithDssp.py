import pandas as pd
import seaborn as sns
import DsspTrainer as dssp

fraction = 0.5
recreate = True

csvfiles = ["Csv/StructureCompareA.csv","Csv/StructureCompareB.csv"]

for csvfile in csvfiles:
    print('### training classifier ########################')
    trainingset = pd.read_csv('Csv/SSTrainingSet_hb.csv')
    trainer = dssp.DsspTrainer(trainingset,fraction=fraction,geos=['N:CA:C:N+1','C-1:N:CA:C'])

    dssp_csv = csvfile + '_dssp.csv'
    print('reading', csvfile)
    data = pd.read_csv(csvfile)

    print('### Creating dataframe ########################')
    data_dssp = trainer.classifySecondaryStructures(data)
    print(data_dssp)
    print('Save to', dssp_csv)
    data_dssp.to_csv(dssp_csv, index=False)





