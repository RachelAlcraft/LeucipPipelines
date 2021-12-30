from os import listdir
from os.path import isfile, join
import os

def unzipProteome(directory):
    print('Unzipping ',directory)
    dir = directory
    onlyfiles = [f for f in listdir(dir) if isfile(join(dir, f))]

    pdbdel = []
    for i in range(0, len(onlyfiles)):
        fn = onlyfiles[i]
        print('nrm',i, '/', len(onlyfiles))
        import gzip
        import shutil
        with gzip.open(dir + fn, 'rb') as f_in:
            fndot = fn.split('.')
            fin = dir + fn
            fout = dir + fndot[0] + '.pdb'
            if fndot[1] == 'cif':
                pdbdel.append(fin)
                print('Add del:', fin)
            elif len(fndot) > 2:
                if fndot[2] == 'gz':
                    pdbdel.append(fin)
                    print('Add del:', fin)
                    try:
                        with open(fout, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    except:
                        print('Error unzipping del:', fin)




    # 2) Delete zipped ###
    print('2) Deleting ###############################')
    for fd in pdbdel:
        os.remove(fd)
