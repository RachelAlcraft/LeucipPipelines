


import numpy as np
from LeucipPy import GeoHTML as ghm
from LeucipPy import LeuciFile as lef

def extractPdb(LogFile, pdbcodes):
    for pdbcode in pdbcodes:
        lfil = lef.LeuciFile(LogFile + pdbcode + "_EMBELLISHFILE.csv")
        new_file = LogFile + "pdb_embellish" + pdbcode + ".ent"
        lfil.saveText('EMBELLISH',new_file)