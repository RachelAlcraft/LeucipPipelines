
run,extract = True,True
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
pdbcodes = ['1aho']
ExePath ='C:/Dev/Github/PsuMaxima/Linux/out/build/x64-Release/PsuMaxima.exe'
Ccp4Directory = "C:/Dev/Github/ProteinDataFiles/ccp4_data/"
PdbDirectory = "C:/Dev/Github/ProteinDataFiles/pdb_data/"
LogDirectory = "C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Log/"
LogDirectory = "Csv/"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import A01RunLeucipPlus as A01
import A02ExtractPdb as A02
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if run:
    A01.runCreateElectronDensity(pdbcodes,ExePath,LogDirectory)
if extract:
    A02.extractPdb(LogDirectory,pdbcodes)