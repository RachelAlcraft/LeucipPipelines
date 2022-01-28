
import os
import subprocess as sub
import shutil


def runGemmi(pdbs, dataDir, pdbDir, ccp4Dir):
    '''
    https://gemmi.readthedocs.io/en/latest/utils.html#sf2map
    :return:
    '''
    '''
    :return: 
    '''
    #import gemmi
    #mtz_g = gemmi.read_mtz_file(mtz_fullfile)
    #gemmi.Mtz.transform_f_phi_to_map()
    # or use the exe
    gemmiexe = "C:/CCP4-7/7.1/bin/gemmi.exe"
    cmda = "sf2map"

    for pdb in pdbs:
        print(pdb,'redo mtz conversion')
        pdbf = dataDir + "pdb" + pdb + ".ent"
        cmdb = dataDir + "pdb" + pdb + ".mtz"
        cmdc = dataDir + pdb + ".ccp4"
        pigP = sub.Popen([gemmiexe, cmda,cmdb,cmdc], stdout=sub.PIPE)
        resultP = pigP.communicate(input=b"This is sample text.\n")
        exe_resultP = str(resultP[0], 'utf-8')
        pigP.kill()

        shutil.copy(pdbf, pdbDir + "/pdb" + pdb + ".ent")
        shutil.copy(cmdc, ccp4Dir + "/" + pdb + ".ccp4")
        shutil.copy(cmdc, ccp4Dir + "/" + pdb + "_diff.ccp4")
