'''
https://gemmi.readthedocs.io/en/latest/utils.html#sf2map
'''
import subprocess as sub
import shutil

pdb = '1aho'

#inuts from pdb or pdb redo
orig_ccp4 = 'C:/Dev/Github/ProteinDataFiles/ccp4_data/' + pdb + '.ccp4'
orig_mtz = 'C:/Dev/Github/ProteinDataFiles/mtz_data/' + pdb + '_phases.mtz'
redo_mtz =  'C:/Dev/Github/ProteinDataFiles/mtz_data/' + pdb + '_final.mtz'

#gemmi conversions
gemmi_dir = 'C:/Dev/Github/ProteinDataFiles/gemmi_conversions/'
orig_converted_to_mtz = gemmi_dir + pdb + '_orig_ccp4.mtz'
orig_back_to_ccp4 = gemmi_dir + pdb + '_orig_mtz.ccp4'
orig_converted_to_mmcif = gemmi_dir + pdb + '_orig_ccp4.mmcif'
orig_converted_mtz_to_mmcif = gemmi_dir + pdb + '_orig_mtz.mmcif'
redo_converted_to_mmcif = gemmi_dir + pdb + '_redo_mtz.mmcif'
redo_converted_to_ccp4 = gemmi_dir + pdb + '_redo_mtz.ccp4'
redo_converted_ccp4_mtz = gemmi_dir + pdb + '_redo_ccp4.mtz'


gemmiexe = "C:/CCP4-7/7.1/bin/gemmi.exe"

commandsA = []
commandsA.append(['map2sf',orig_ccp4,orig_converted_to_mtz,'FWT','PHWT'])

commandsB = []
commandsB.append(['mtz2cif',orig_mtz,orig_converted_mtz_to_mmcif])
commandsB.append(['mtz2cif',redo_mtz,redo_converted_to_mmcif])
#reconversions
commandsB.append(['mtz2cif',orig_converted_to_mtz,orig_converted_to_mmcif])

commandsC = []
commandsC.append(['sf2map',orig_converted_to_mtz,orig_back_to_ccp4,'--grid=-47,-54,-62'])

for cmd,a,b,c,d in commandsA:
    pigP = sub.Popen([gemmiexe, cmd,a,b,c,d], stdout=sub.PIPE)
    resultP = pigP.communicate(input=b"This is sample text.\n")
    exe_resultP = str(resultP[0], 'utf-8')
    pigP.kill()

for cmd,a,b in commandsB:
    pigP = sub.Popen([gemmiexe, cmd,a,b], stdout=sub.PIPE)
    resultP = pigP.communicate(input=b"This is sample text.\n")
    exe_resultP = str(resultP[0], 'utf-8')
    pigP.kill()

for cmd,a,b,c in commandsC:
    pigP = sub.Popen([gemmiexe, cmd,a,b,c], stdout=sub.PIPE)
    resultP = pigP.communicate(input=b"This is sample text.\n")
    exe_resultP = str(resultP[0], 'utf-8')
    pigP.kill()


