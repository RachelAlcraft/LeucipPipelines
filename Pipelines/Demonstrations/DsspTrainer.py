

'''
Training for DSSP
'''
import pandas as pd
from sklearn.model_selection import cross_validate
from sklearn.ensemble import RandomForestClassifier


class DsspTrainer:
    def __init__(self,trainingset,fraction=1,log=0,
                 geos = ['N:O+2', 'N:O+3', 'N:O+4', 'O:N+2', 'O:N+3', 'O:N+4', 'O:{N}+2', 'N:{O}+2', 'N:CA:C', 'N:CA:C:N+1','C-1:N:CA:C'],
                 aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN','ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']):
        '''
        Load up the training set and train in preparation
        A useful link to dssp: https://ssbio.readthedocs.io/en/latest/instructions/dssp.html
        H	Alpha helix
        B	Beta bridge
        E	Strand
        G	Helix-3
        I	Helix-5
        T	Turn
        S	Bend
        H,0.9140212099475761
        T,0.8260132320102317
        E,0.8256354104302333
        S,0.5476708182191182
        G,0.5387212643678161
        B,0.5
        I,0.5

        First geos were:
        geos=['N:CA','CA:C','C:N+1','C:O','N:CA:C:N+1','C-1:N:CA:C','CA-1:C-1:N:CA','CA:C:N+1:CA+1','CA-1:CA:CA+1','CA-1:CA+1','CA:CA+1','CA:CA-1'],
        '''
        pd.options.mode.chained_assignment = None
        self.fraction = fraction
        self.trainingset = trainingset#pd.read_csv('Csv/dssp_trainingset.csv')

        if log > 0:
            datagrouped = self.trainingset[['dssp', 'aa']]
            datagrouped = datagrouped.groupby(by=['dssp']).count()
            datagrouped['dssp'] = datagrouped.index
            datagrouped['count'] = datagrouped['aa']
            print('LeucipPy (1) dssp split=',datagrouped[['dssp', 'count']])


        print(self.trainingset)
        #I want to have a training set for every amino acid AND for every dssp
        self.aa_list = aa_list
        #self.aa_list = ['ALL']
        #self.aa_list = ['GLY']
        self.dssp_list = ['H','T','E','S','G','B','I']
        self.dssp_tag = {'H':'AlphaHelix','T':'Turn','E':'Strand','S':'Bend','G':'Helix3','B':'BetaBridge'
                                                                                              '','I':'Helix5'}
        self.geos = geos
        self.trainers = {}
        for aa in self.aa_list:
            self.trainers[aa] = {}
            for dssp in self.dssp_list:
                print(aa, dssp,'...',end='')
                if aa != 'ALL':
                    dataaa = self.trainingset.query("aa == '" + aa + "'").sample(frac=fraction)
                else:
                    dataaa = self.trainingset.sample(frac=fraction)
                sss = []
                for i in range(len(dataaa.index)):
                    row = dataaa.iloc[i,:]
                    if row.loc['dssp'] == dssp:
                        sss.append(1)
                    else:
                        sss.append(0)
                dataaa.loc[:,'ss'] = sss
                dataaa = dataaa.drop('dssp',axis=1)
                dataaa = dataaa.drop('aa', axis=1)

                # manipulate data
                y = dataaa.loc[:,'ss'].values
                X = dataaa[self.geos].values

                # train classifier
                metrics = ['balanced_accuracy']
                clf = RandomForestClassifier(random_state=1)
                s = cross_validate(clf, X, y, cv=2, scoring=metrics, return_train_score=False)
                print('...',s['test_balanced_accuracy'].mean(),end='\n')
                clf.fit(X,y)
                self.trainers[aa][dssp] = clf

    def classifySecondaryStructures(self,data):
        '''
        :param data: it must have the geos used for secondary structure prediction
        :return: a new dataframe with an extra column ss
        '''
        predictions = []
        for i in range(len(data.index)):
            row = data.iloc[i, :]
            rowaa = row['aa']
            if rowaa in  self.aa_list:
                if i%100 == 0:
                    print(i,'/',len(data.index),rowaa)
                rowX = row[self.geos]
                classifiers = self.trainers[rowaa]
                #classifiers = self.trainers['ALL']
                found = False
                for ss,clf in classifiers.items():
                    #print(ss,clf)
                    if not found:
                        predict = clf.predict([rowX.values])
                        #print('Prediction=',predict)
                        if predict == 1:
                            predictions.append(self.dssp_tag[ss])
                            found = True
                if not found:
                    predictions.append('-')

            else:
                predictions.append('-')

        #print(predictions)
        data['ss'] = predictions
        return data


