'''
Generating positions with external library
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3661355/
'''

from LeucipPy import GeometryMaker as gm
from LeucipPy import HtmlReportMaker as hm
import pandas as pd
import Ext_Geometry as ext_geo
import Ext_PeptideBuilder as ext_pep
import A0Class_AtomCollection as rot

omega_rules = rot.RotationRules('0.5{-180,-150}:0.5{150,180}')
phi_rules = rot.RotationRules('0.5{-180,-50}:0.5{50,180}')
psi_rules = rot.RotationRules('4{-180,-120}:1{-120,-50}:4{-50,50}:1{50,120}:4{120,180}')

strucs = []

for i in range(0,1000):
    omega = int(omega_rules.getRandomRotation())
    phi = int(phi_rules.getRandomRotation())
    psi = int(psi_rules.getRandomRotation())
    print(omega,phi,psi)
    geo = ext_geo.geometry('G')
    geo.omega = omega
    geo.phi = phi
    geo.psi_im1 = psi
    structure = ext_pep.initialize_res(geo)
    structure.header = {}
    structure.header['resolution'] = 1
    structure.id = 'X' + str(i)

    for i in range(0,2):
        structure = ext_pep.add_residue(structure,geo)
    strucs.append(structure)

geo_mak = gm.GeometryMaker(strucs)

df = geo_mak.calculateGeometry(['N:CA:C:N+1','N:O','C-1:N:CA:C','N:CA','CA:C','C:O','C:N+1','CA:C:N+1:CA+1','CA-1:C-1:N:CA','N:N+1','C-1:C','O:CA-1','O:C-1','O-1:N+1','C-1:N+1','CA-1:CA:CA+1','CA-1:CA+1'])
print(df)
df.to_csv('Csv/A0_SynthBio.csv',index=False)
rep_mak = hm.HtmlReportMaker('Synthetic Geometry','Html/AO_SyntheticRandomBioMaker.html', cols=4)
rep_mak.addPlot1d(df, 'histogram', 'C-1:N:CA:C', hue='pdb_code', title='PHI')
rep_mak.addPlot1d(df, 'histogram', 'N:CA:C:N+1', hue='pdb_code', title='PSI')
rep_mak.addPlot1d(df, 'histogram', 'CA:C:N+1:CA+1', hue='pdb_code', title='OMEGA')
rep_mak.addPlot1d(df, 'histogram', 'N:N+1', hue='pdb_code', title='N:N')

rep_mak.addPlot2d(df, 'scatter', 'C-1:N:CA:C', 'N:CA:C:N+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'C-1:N:CA:C', 'C-1:C', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'C-1:N+1', 'O-1:N+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'CA-1:CA+1', 'CA-1:CA:CA+1', hue='N:CA:C:N+1', title='', palette='jet_r')


rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1', 'C-1:N:CA:C', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1', 'CA:C:N+1:CA+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1', 'CA-1:C-1:N:CA', hue='N:CA:C:N+1', title='', palette='jet_r')

rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1','N:O', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.addPlot2d(df, 'scatter', 'N:CA:C:N+1','N:N+1', hue='N:CA:C:N+1', title='', palette='jet_r')
rep_mak.printReport()