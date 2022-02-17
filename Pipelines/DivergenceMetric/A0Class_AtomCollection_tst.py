
import A0Class_AtomCollection as syn

# make the atoms dataframe to pass in
import pandas as pd

dic_atoms = {}
dic_atoms['branch'] = ['1','1','1','1.2']
dic_atoms['position'] =['1','2','3','4.1']
dic_atoms['resNo'] =[1,1,1,1]
dic_atoms['aminoAcid'] =['GLY','GLY','GLY','GLY']
dic_atoms['atom']=['N','CA','C','O']
dic_atoms['coordX'] =[0,1,2,3]
dic_atoms['coordY'] =[0,0,0,0]
dic_atoms['coordZ'] =[0,0,0,0]
dic_atoms['rotationSet'] =['1{0,90}','1{0,90}','1{0,90}','{}']
atoms_data = pd.DataFrame.from_dict(dic_atoms)
atoms_maker = syn.AtomCollection(atoms_data,log=2)
print(atoms_maker)
atoms_df = atoms_maker.makeRandomVersions(2)
print(atoms_df)


