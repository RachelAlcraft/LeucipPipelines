'''
RSA 28-Feb-22
'''
import math
import random
import Ext_Geometry as ext_geo
import numpy as np
import pandas as pd
from LeucipPy import GeoCalculate as calc

class AtomCoordinates:
    def __init__(self,data, log=0):
        self.data = data
        self.log = log
        self.rand =random.randint(0,len(self.data.index)-1)
        if log > 1:
            print('LeucipPy(2), rand=',self.rand)

    def getRandomParams(self,geos,regen=True):
        ret_dic = {}
        if regen:
            self.rand =random.randint(0,len(self.data.index)-1)
            if self.log > 1:
                print('LeucipPy(2), rand=', self.rand)
        for geo in geos:
            val = self.data[geo].values[self.rand]
            ret_dic[geo] = float(val)

        return ret_dic

    def generateAllRandomGeo(self):
        geom = ext_geo.geometry('G')
        # dihedrals
        geom.omega = self.getRandomParams(['CA-1:C-1:N:CA'],regen=True)['CA-1:C-1:N:CA']
        geom.phi = self.getRandomParams(['C-1:N:CA:C'],regen=True)['C-1:N:CA:C']
        geom.psi_im1 = self.getRandomParams(['N-1:CA-1:C-1:N'],regen=True)['N-1:CA-1:C-1:N']
        psi_next = self.getRandomParams(['N:CA:C:N+1'], regen=False)['N:CA:C:N+1']
        geom.N_CA_C_O_diangle = self.getRandomParams(['N:CA:C:O'],regen=True)['N:CA:C:O']
        # angles =
        geom.N_CA_C_angle = self.getRandomParams(['N:CA:C'],regen=True)['N:CA:C']
        geom.C_N_CA_angle = self.getRandomParams(['C-1:N:CA'],regen=True)['C-1:N:CA']
        geom.CA_C_N_angle = self.getRandomParams(['CA-1:C-1:N'],regen=True)['CA-1:C-1:N']
        geom.CA_C_O_angle = self.getRandomParams(['CA:C:O'],regen=True)['CA:C:O']
        # bonds =
        geom.CA_N_length = self.getRandomParams(['N:CA'],regen=True)['N:CA']
        geom.CA_C_length = self.getRandomParams(['CA:C'],regen=True)['CA:C']
        geom.C_O_length = self.getRandomParams(['C:O'],regen=True)['C:O']
        geom.peptide_bond = self.getRandomParams(['N:C-1'],regen=True)['N:C-1']
        return geom, psi_next

    def generateSampleGeo(self):
        geom = ext_geo.geometry('G')
        # dihedrals
        geom.omega = self.getRandomParams(['CA-1:C-1:N:CA'],regen=False)['CA-1:C-1:N:CA']
        geom.phi = self.getRandomParams(['C-1:N:CA:C'],regen=False)['C-1:N:CA:C']
        geom.psi_im1 = self.getRandomParams(['N-1:CA-1:C-1:N'],regen=False)['N-1:CA-1:C-1:N']
        psi_next = self.getRandomParams(['N:CA:C:N+1'], regen=False)['N:CA:C:N+1']
        geom.N_CA_C_O_diangle = self.getRandomParams(['N:CA:C:O'],regen=False)['N:CA:C:O']
        # angles =
        geom.N_CA_C_angle = self.getRandomParams(['N:CA:C'],regen=False)['N:CA:C']
        geom.C_N_CA_angle = self.getRandomParams(['C-1:N:CA'],regen=False)['C-1:N:CA']
        geom.CA_C_N_angle = self.getRandomParams(['CA-1:C-1:N'],regen=False)['CA-1:C-1:N']
        geom.CA_C_O_angle = self.getRandomParams(['CA:C:O'],regen=False)['CA:C:O']
        # bonds =
        geom.CA_N_length = self.getRandomParams(['N:CA'],regen=False)['N:CA']
        geom.CA_C_length = self.getRandomParams(['CA:C'],regen=False)['CA:C']
        geom.C_O_length = self.getRandomParams(['C:O'],regen=False)['C:O']
        geom.peptide_bond = self.getRandomParams(['N:C-1'],regen=False)['N:C-1']

        return geom,psi_next

    def generateRamaGeo(self):
        geom = ext_geo.geometry('G')
        # dihedrals
        geom.phi = self.getRandomParams(['C-1:N:CA:C'], regen=True)['C-1:N:CA:C']
        geom.psi_im1 = self.getRandomParams(['N-1:CA-1:C-1:N'], regen=False)['N-1:CA-1:C-1:N']
        psi_next = self.getRandomParams(['N:CA:C:N+1'], regen=False)['N:CA:C:N+1']
        geom.omega = self.getRandomParams(['CA-1:C-1:N:CA'],regen=True)['CA-1:C-1:N:CA']
        geom.N_CA_C_O_diangle = self.getRandomParams(['N:CA:C:O'],regen=True)['N:CA:C:O']
        # angles =
        geom.N_CA_C_angle = self.getRandomParams(['N:CA:C'],regen=True)['N:CA:C']
        geom.C_N_CA_angle = self.getRandomParams(['C-1:N:CA'],regen=True)['C-1:N:CA']
        geom.CA_C_N_angle = self.getRandomParams(['CA-1:C-1:N'],regen=True)['CA-1:C-1:N']
        geom.CA_C_O_angle = self.getRandomParams(['CA:C:O'],regen=True)['CA:C:O']
        # bonds =
        geom.CA_N_length = self.getRandomParams(['N:CA'],regen=True)['N:CA']
        geom.CA_C_length = self.getRandomParams(['CA:C'],regen=True)['CA:C']
        geom.C_O_length = self.getRandomParams(['C:O'],regen=True)['C:O']
        geom.peptide_bond = self.getRandomParams(['N:C-1'],regen=True)['N:C-1']
        #geom.peptide_bond = self.getRandomParams(['C:N+1'], regen=True)['C:N+1']
        return geom,psi_next




