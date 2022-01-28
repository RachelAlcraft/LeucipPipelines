
pdbName = '7a4m'
#pdbName = '6eex'
emName = 'EMD-11638'
dir = 'C:/Dev/Github/LeucipPipelines/Pipelines/DeformationDensity/1EMData/Data/'
ExePath ='C:/Dev/Github/PsuMaxima/Linux/out/build/x64-Release/PsuMaxima.exe'


import A01DownloadFiles as A01
import A02CreateTheoreticalIAM as A02

A01.downloadEMFiles(pdbName,emName)
A02.runCreateIAM(pdbName,ExePath,dir)