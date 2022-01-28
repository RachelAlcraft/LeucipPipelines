import gzip
import os
import shutil
from urllib.request import urlretrieve


def downloadEMFiles(pdbName, emName):
    em_new = emName.replace('-','_').lower()
    pdb_path = 'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdbName + '.ent'
    em_path = 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/' +emName+ '/map/' +em_new+ '.map.gz'
    pdb_file = 'C:/Dev/Github/LeucipPipelines/Pipelines/DeformationDensity/1EMData/Data/pdb' + pdbName + '.ent'
    em_file_gz = 'C:/Dev/Github/LeucipPipelines/Pipelines/DeformationDensity/1EMData/Data/' + pdbName + '.ccp4.gz'
    em_file = 'C:/Dev/Github/LeucipPipelines/Pipelines/DeformationDensity/1EMData/Data/' + pdbName + '.ccp4'
    em_file_diff = 'C:/Dev/Github/LeucipPipelines/Pipelines/DeformationDensity/1EMData/Data/' + pdbName + '_diff.ccp4'

    if not os.path.exists(pdb_file):
        try:
            urlretrieve(pdb_path, pdb_file)
        except:
            print("...No pdb data for", pdbName)

    if not os.path.exists(em_file):
        try:
            print(em_path,em_file_gz)
            urlretrieve(em_path, em_file_gz)
            with gzip.open(em_file_gz, 'rb') as f_in:
                with open(em_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                with open(em_file_diff, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        except:
            print("...No map data for", pdbName)

