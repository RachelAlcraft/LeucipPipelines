
import A0Class_AtomCollection as syn
from LeucipPy import GeometryMaker as gm
from LeucipPy import HtmlReportMaker as hm
# make the atoms dataframe to pass in
import pandas as pd

#taken 11,12,13 N,CA,C,O coords from 1aba
'''
						x	y	z
ATOM	177	N			1	7.224	2.499	8.39
ATOM	178	CA			1	6.524	1.546	7.525
ATOM	179	C			1	5.418	2.215	6.717
ATOM	180	O			1	5.249	1.974	5.505
ATOM	187	N			2	4.655	3.072	7.407
ATOM	188	CA			2	3.604	3.824	6.745
ATOM	189	C			2	4.15	4.797	5.7
ATOM	190	O			2	3.537	4.954	4.628
ATOM	201	N			3	5.248	5.433	6.024
ATOM	202	CA			3	5.942	6.331	5.086
ATOM	203	C			3	6.258	5.59	3.779
ATOM	204	O			3	5.99	6.074	2.638
'''

dic_atoms = {'branch':[],'position':[],'resNo':[],'aminoAcid':[],'atom':[],'coordX':[],'coordY':[],'coordZ':[],'rotationSet':[]}
#rotation_set = ['{}','1{200,210}:1{100,110}:1{350,360}','{}','{}'] #PHI
rotation_set = ['{}','1{10,10}',  '{}','{}'] #PHI
#rotation_set = ['{}','{}','1{0,360}','{}'] #PSI
#rotation_set = ['{}','1{0,360}','1{0,360}','{}'] #both
#rotation_set = ['{}','{}','{}','1{0,360}'] #CO!
#rotation_set = ['0.5{0,20}:0.5{340,360}','{}','{}','{}'] #omega
#rotation_set = ['0.5{0,20}:0.5{340,360}','1{0,360}','1{0,360}','{}'] #ALL
#rotation_set = ['{}','{}','{}','{}']
# Residue 1
dic_atoms['branch'].extend(['1','1','1','1.2'])
dic_atoms['position'].extend(['1','2','3','3.1'])
dic_atoms['resNo'].extend([1,1,1,1])
dic_atoms['aminoAcid'].extend(['GLY','GLY','GLY','GLY'])
dic_atoms['atom'].extend(['N','CA','C','O'])
dic_atoms['coordX'].extend([7.224,6.524,5.418,5.249])
dic_atoms['coordY'].extend([2.499,1.546,2.215,1.974])
dic_atoms['coordZ'].extend([8.39,7.525,6.717,5.505])
#dic_atoms['rotationSet'].extend(['0.5{0,40}:0.5{320,360}','1{0,360}','1{0,360}','{}'])
dic_atoms['rotationSet'].extend(rotation_set)
# Residue 2
dic_atoms['branch'].extend(['1','1','1','1.3'])
dic_atoms['position'].extend(['4','5','6','6.1'])
dic_atoms['resNo'].extend([2,2,2,2])
dic_atoms['aminoAcid'].extend(['GLY','GLY','GLY','GLY'])
dic_atoms['atom'].extend(['N','CA','C','O'])
dic_atoms['coordX'].extend([4.655,3.604,4.15,3.537])
dic_atoms['coordY'].extend([3.072,3.824,4.797,4.954])
dic_atoms['coordZ'].extend([7.407,6.745,5.7,4.628])
#dic_atoms['rotationSet'].extend(['0.5{0,40}:0.5{320,360}','1{0,360}','1{0,360}','{}'])
dic_atoms['rotationSet'].extend(rotation_set)
# Residue 3
dic_atoms['branch'].extend(['1','1','1','1.4'])
dic_atoms['position'].extend(['7','8','9','9.1'])
dic_atoms['resNo'].extend([3,3,3,3])
dic_atoms['aminoAcid'].extend(['GLY','GLY','GLY','GLY'])
dic_atoms['atom'].extend(['N','CA','C','O'])
dic_atoms['coordX'].extend([5.248,5.942,6.258,5.99])
dic_atoms['coordY'].extend([5.433,6.331,5.59,6.074])
dic_atoms['coordZ'].extend([6.024,5.086,3.779,2.638])
#dic_atoms['rotationSet'].extend(['0.5{0,40}:0.5{320,360}','1{0,360}','1{0,360}','{}'])
dic_atoms['rotationSet'].extend(rotation_set)

atoms_data = pd.DataFrame.from_dict(dic_atoms)
atoms_maker = syn.AtomCollection(atoms_data,log=2)
print(atoms_maker)
atoms_df = atoms_maker.makeRandomVersions(10)
print(atoms_df)

geo_mak = gm.GeometryMaker([],init_biopython=False,atoms_data=atoms_df)
df = geo_mak.calculateGeometry(['N:CA:C:N+1','N:O','C-1:N:CA:C','N:CA','CA:C','C:O','C:N+1','CA:C:N+1:CA+1','CA-1:C-1:N:CA','N:N+1','C-1:C'])
print(df)

rep_mak = hm.HtmlReportMaker('Synthetic Geometry','Html/AO_SyntheticRandom.html', cols=5)
# Some scatter plots are added to make better sense of the data, the colour scheme is to approcimate the AlphaFold probabilty
#rep_mak.addPlot1d(df, 'histogram', 'N:CA', hue='pdb_code', title='')
#rep_mak.addPlot1d(df, 'histogram', 'CA:C', hue='pdb_code', title='')
#rep_mak.addPlot1d(df, 'histogram', 'C:O', hue='pdb_code', title='')
#rep_mak.addPlot1d(df, 'histogram', 'C:N+1', hue='pdb_code', title='')
rep_mak.addPlot1d(df, 'histogram', 'C-1:N:CA:C', hue='pdb_code', title='PHI')
rep_mak.addPlot1d(df, 'histogram', 'N:CA:C:N+1', hue='pdb_code', title='PSI')
rep_mak.addPlot1d(df, 'histogram', 'CA:C:N+1:CA+1', hue='pdb_code', title='OMEGA')
rep_mak.addPlot1d(df, 'histogram', 'N:N+1', hue='pdb_code', title='N:N')
rep_mak.addPlot1d(df, 'histogram', 'C-1:C', hue='pdb_code', title='C:C')
rep_mak.addPlot1d(df, 'histogram', 'N:O', hue='pdb_code', title='N:O')
rep_mak.addPlot1d(df, 'histogram', 'N:CA', hue='pdb_code', title='N:CA')
rep_mak.addPlot1d(df, 'histogram', 'CA:C', hue='pdb_code', title='CA:C')
rep_mak.addPlot1d(df, 'histogram', 'C:O', hue='pdb_code', title='C:O')
#rep_mak.addPlot2d(df, 'scatter', 'N:CA', 'CA:C', hue='N:CA', title='', palette='jet_r')
#rep_mak.addPlot2d(df, 'scatter', 'C:N+1', 'C:O', hue='C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'C-1:N:CA:C', 'CA:C:N+1:CA+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'C-1:N:CA:C', 'CA-1:C-1:N:CA', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1','N:O', hue='N:CA:C:N+1', title='', palette='jet_r')

rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1', 'C-1:N:CA:C', hue='C-1:N:CA:C', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1', 'CA:C:N+1:CA+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1', 'CA-1:C-1:N:CA', hue='N:CA:C:N+1', title='', palette='jet_r')

rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1','N:N+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.printReport()



