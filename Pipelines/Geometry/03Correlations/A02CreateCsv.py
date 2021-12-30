from os import listdir
from os.path import isfile, join
import os
from LeucipPy import GeoDataFrame as gdf
import Bio.PDB as bio
import pandas as pd


def createCsvCorrelations(name,directory,geos,ext,chunk,rewrite):
    print('Creating correlations from ',directory)
    dir = directory
    csv_final = "Csv/" + name + "_" + ext + ".csv"

    if not os.path.exists(csv_final) or rewrite:
        onlyfiles = [f for f in listdir(dir) if isfile(join(dir, f))]
        parser = bio.PDBParser()
        dfs = []
        for c in range(0, len(onlyfiles),chunk):
            start = c
            end = c + chunk
            if end > len(onlyfiles):
                end = len(onlyfiles)
            csv_name = "Csv/" + name + "_" + str(start+1) + '_' + str(end) + "_" + ext + ".csv"
            if not os.path.exists(csv_name) or rewrite:
                print(name,' chunk',start,end)
                strucs = []
                for i in range(start, end):
                    fn = onlyfiles[i]
                    fn = directory + fn
                    print(name, c,':',i, '/', len(onlyfiles),'loading pdb for correlations')
                    fndot = fn.split('.')
                    fnfile = fndot[0].split('/')
                    fnnam = fnfile[len(fnfile) - 1:][0]
                    struc = parser.get_structure(fnnam, fn)
                    strucs.append(struc)
                print(name, c, '/', len(onlyfiles), 'creating dataframe for correlations')
                geo = gdf.GeoDataFrame(strucs, log=0)
                data = geo.calculateGeometry(geos, log=0)
                data.to_csv(csv_name, index=False)
                dfs.append(data)
            else:
                data = pd.read_csv(csv_name)
                dfs.append(data)

        # 4) Create a dataframe of geos
        print('Concatenating dataframe',name)
        df_geometry = pd.concat(dfs)

        print('5) Save dataframe ###############################')
        df_geometry.to_csv(csv_final , index=False)
        print("Saved to", csv_final)


def createCsvCorrelationsFromList(name,directory,geos,ext,pdbs,chunk,rewrite):
    print('Creating correlations from ',directory)
    csv_final = "Csv/" + name + "_" + ext + ".csv"
    print(csv_final)

    if not os.path.exists(csv_final) or rewrite:
        parser = bio.PDBParser()
        dfs = []
        for c in range(0, len(pdbs),chunk):
            start = c
            end = c + chunk
            if end > len(pdbs):
                end = len(pdbs)
            csv_name = "Csv/" + name + "_" + str(start+1) + '_' + str(end) + "_" + ext + ".csv"
            if not os.path.exists(csv_name) or rewrite:
                print(name,' chunk',start,end)
                strucs = []
                for i in range(start, end):
                    fn = 'pdb' + pdbs[i] + '.ent'
                    fn = directory + fn
                    print(name, c,':',i, '/', len(pdbs),'loading pdb for correlations')
                    fndot = fn.split('.')
                    fnfile = fndot[0].split('/')
                    fnnam = fnfile[len(fnfile) - 1:][0]
                    struc = parser.get_structure(fnnam, fn)
                    strucs.append(struc)
                print(name, c, '/', len(pdbs), 'creating dataframe for correlations')
                geo = gdf.GeoDataFrame(strucs, log=0)
                data = geo.calculateGeometry(geos, log=0)
                data.to_csv(csv_name, index=False)
                dfs.append(data)
            else:
                data = pd.read_csv(csv_name)
                dfs.append(data)

        # 4) Create a dataframe of geos
        print('Concatenating dataframe',name)
        df_geometry = pd.concat(dfs)

        print('5) Save dataframe ###############################')
        df_geometry.to_csv(csv_final , index=False)
        print("Saved to", csv_final)



