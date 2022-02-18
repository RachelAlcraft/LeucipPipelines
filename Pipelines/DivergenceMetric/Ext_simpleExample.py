#!/usr/bin/python

'''
Simple example script demonstrating how to use the PeptideBuilder library.

The script generates a peptide consisting of six arginines in alpha-helix
conformation, and it stores the peptide under the name "example.pdb".
'''


import Ext_Geometry
import Ext_PeptideBuilder


geo = Ext_Geometry.geometry('G')
geo.phi=-60
geo.psi_im1=-40
structure = Ext_PeptideBuilder.initialize_res(geo)
for i in range(5):
    structure = Ext_PeptideBuilder.add_residue(structure, geo)
    
import Bio.PDB
out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save( "example.pdb" )


#My added librarary
from LeucipPy import GeometryMaker as gm

structure.header = {}
structure.header['resolution'] = 1
geo_mak = gm.GeometryMaker([structure])
df = geo_mak.calculateGeometry(['N:CA:C:N+1'])
print(df)
df = geo_mak.calculateGeometry(['N:O','N:CA','CA:C','C:O','C:N+1','N:N+1','C-1:C'])
print(df)
df = geo_mak.calculateGeometry(['C-1:N:CA:C'])
#df = geo_mak.calculateGeometry(['N:CA:C:N+1','N:O','C-1:N:CA:C','N:CA','CA:C','C:O','C:N+1','CA:C:N+1:CA+1','CA-1:C-1:N:CA','N:N+1','C-1:C'])
print(df)
df = geo_mak.calculateGeometry(['CA:C:N+1:CA+1'])
print(df)
