'''
https://gemmi.readthedocs.io/en/latest/utils.html#sf2map
'''
import subprocess as sub
import shutil

### This file transforms an mtz to a cif so we can look inside it ###

#inuts from pdb or pdb redo
mtz = 'C:/Users/Rachel/Google Drive/ProjectPlan/Education/PhDBio/Practicals/001_Tracey01/mug_dna.mtz'
cif = 'C:/Users/Rachel/Google Drive/ProjectPlan/Education/PhDBio/Practicals/001_Tracey01/gemmi_mug_dna_mtz.cif'
pycif = 'C:/Users/Rachel/Google Drive/ProjectPlan/Education/PhDBio/Practicals/001_Tracey01/py_gemmi_mug_dna_mtz.cif'
#gemmi conversions
gemmiexe = "C:/CCP4-7/7.1/bin/gemmi.exe"


commands = [['mtz2cif',mtz,cif]]

for cmd,a,b in commands:
    pigP = sub.Popen([gemmiexe, cmd,a,b], stdout=sub.PIPE)
    resultP = pigP.communicate(input=b"This is sample text.\n")
    exe_resultP = str(resultP[0], 'utf-8')
    pigP.kill()

import sys
from matplotlib import pyplot
import gemmi

mtzobj = gemmi.read_mtz_file(mtz)
print(mtzobj.datasets)
mincol = 10000000
f = open(pycif, "w")
f.write('loop_\n')
line = ""
for col in mtzobj.columns:
    mincol = min(len(col),mincol)
    line += '_' + str(col).split(" ")[1] + "\n"
    print(col)

for i in range(0,mincol):
    f.write(line[:-1] + '\n')
    line = ""
    for col in mtzobj.columns:
        col_val = col[i]
        try:
            if int(col_val) != col_val:
                col_val ="{:.6f}".format(col[i])
            else:
                col_val = int(col_val)
            line += str(col_val) + " "
        except:
            line += str(col_val) + " "

f.write(line[:-1] + '\n')
f.close()
