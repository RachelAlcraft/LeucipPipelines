'''
This manages the geoemtric movement of atoms based on the following dataframe folamr
Branch	Position	ResNo	AminoAcid	Atom	Coords	Rotation
1	1	1	GLY	N	(0,0,0)	{}
1	2	1	GLY	CA	etc	{0,360}
1	3	1	GLY	C	etc	{0,360}
1.2.1	4.1	1	GLY	O	etc	{}
1	4	2	GLY	N	etc	{0,40},{320,360}
1	5	2	GLY	CA	etc	{0,360}
1	6	2	GLY	C	etc	{0,360}
1.2.2	7.1	2	GLY	O	etc	{}
1	7	3	GLY	N	etc	{0,40},{320,360}
1	8	3	GLY	CA	etc	{0,360}
1	9	3	GLY	C	etc	{0,360}
1.2.3	10.1	3	GLY	O	etc	{}
'''
import math
import random

import pandas as pd
from LeucipPy import GeoCalculate as calc

M_PI = 3.14159265358979323846   # pi


class RotationRules:
    def __init__(self,rotationString):
        self.prob_list = []
        self.prob_hat = []
        self.range_list = []
        if rotationString != '{}':
            sets = rotationString.split(':')
            for st in sets:
                parts = st.split('{')
                prob = parts[0]
                angle_range = parts[1].split(',')
                rangeA = angle_range[0]
                rangeB = angle_range[1][:-1]
                self.prob_list.append(float(prob))
                self.range_list.append([rangeA,rangeB])

        for i in range(0,len(self.prob_list)):
            prob = self.prob_list[i]
            vals = int(prob*100)
            for v in range(0,vals):
                self.prob_hat.append(i)
        #print(self.prob_hat)

    def getRandomRotation(self):
        # first randomly choose which probaility list to take from
        if len(self.prob_hat) > 0:
            random.shuffle(self.prob_hat)
            hat = self.prob_hat[0]
            rangeA = int(self.range_list[hat][0])
            rangeB = int(self.range_list[hat][1])
            angle_hat = []
            for i in range(rangeA,rangeB+1):
                angle_hat.append(i)
            random.shuffle(angle_hat)
            angle = angle_hat[0]
            return angle
        return None



class Atom:
    def __init__(self, branch,position,resNo,aminoAcid,atom,coordX,coordY,coordZ,rotationSet,mainBranch,link_to_last):
        self.branch,self.position = branch,position
        self.resNo,self.aminoAcid,self.atom, = resNo,aminoAcid,atom
        self.coordX, self.coordY, self.coordZ = coordX,coordY,coordZ
        self.rotation = RotationRules(rotationSet)
        self.mainBranch = mainBranch
        self.last_atom = link_to_last

class AtomCollection:
    def __init__(self,atoms_data,log=0):
        self.atoms_data = atoms_data
        self.log = log
        self.branches_sets = {}


    def moveDataFrameToReadyLists(self):
        # when we init we turn this into sets and a list of atoms, tuples are by ref so all the sets should be by ref
        rids, aas, xs, ys, zs, atoms = [], [], [], [], [], []
        self.branches_sets = {}
        last_atm = None
        for i in range(len(self.atoms_data.index)):
            branches = self.atoms_data['branch'].values[i].split('.')
            positions = self.atoms_data['position'].values[i].split('.')
            resNo = self.atoms_data['resNo'].values[i]
            aminoAcid = self.atoms_data['aminoAcid'].values[i]
            atom = self.atoms_data['atom'].values[i]
            coordX = self.atoms_data['coordX'].values[i]
            coordY = self.atoms_data['coordY'].values[i]
            coordZ = self.atoms_data['coordZ'].values[i]
            rotationSet = self.atoms_data['rotationSet'].values[i]
            mainBranch = len(branches)==1
            atm = Atom(branches[0],positions[0],resNo,aminoAcid,atom,coordX,coordY,coordZ,rotationSet,mainBranch,last_atm)
            if mainBranch:
                last_atm = atm
            for b in range(0,len(branches)):
                if branches[b] not in self.branches_sets:
                    self.branches_sets[branches[b]] = []
                self.branches_sets[branches[b]].append(atm)

    def shiftToPole(self, coordsPole, coordsCentre):
        '''This finds the transformatinios needed to:
        a) Tranlate coords 2 to origin
        b) Rotate coords1 to XZ plane
        c) Rotate coords1 to z axis pole

        a) Translate coords1 to origin
        b) Rotate coords2 to XZ plane
        c) Rotate coords2 to z axis pole
        d) Tranlate coords 2 to origin

        :param coords1: the cords of the first atom
        :param coords2: the coords of the atom that will be rotated
        :return: a tuple of these transforamtions
        '''
        #print('This shifts the coords1 to origin, orients coords2 onto the x axis, then slifers coords2 down')
        transformations = []
        if self.log > 1:
            print('LeucipPy(2) Create Pole Shift for',coordsPole,coordsCentre)
        #a) Tranlate to origin, is just the reverse of the coords
        #a_tranlate = [-1*coords1[0],-1*coords1[1],-1*coords1[2]]
        a_translate = [-1 * coordsCentre[0], -1 * coordsCentre[1], -1 * coordsCentre[2]]
        transformations.append(a_translate)
        coordsNew = [coordsPole[0]+a_translate[0],coordsPole[1]+a_translate[1],coordsPole[2]+a_translate[2]]

        #b) Rotate to y=0 on xy plane
        b_theta_radians = self.__getRotationAngle(coordsNew[0],coordsNew[1])
        b_theta_degrees = self.__radiansToDegrees(b_theta_radians)
        transformations.append([b_theta_degrees])
        coordsNew[0],coordsNew[1] = self.rotate(coordsNew[0],coordsNew[1],b_theta_degrees)

        # c) Rotate to z=0 on XZ plane
        c_theta_radians = self.__getRotationAngle(coordsNew[0], coordsNew[2])
        c_theta_degrees = self.__radiansToDegrees(c_theta_radians)
        transformations.append([c_theta_degrees])
        coordsNew[0], coordsNew[2] = self.rotate(coordsNew[0], coordsNew[2], c_theta_degrees)

        if self.log > 1:
            print('LeucipPy(2) Shifted pole=', round(coordsNew[0],4),round(coordsNew[1],4),round(coordsNew[2],4))

        '''
        # d) Tranlate to origin, is just the reverse of the coords
        d_tranlate = [-1 * coordsNew[0], -1 * coordsNew[1], -1 * coordsNew[2]]
        transformations.append(d_tranlate)        
        '''
        if self.log > 1:
            print('LeucipPy(2) Pole Shift=', transformations)

        # TEST it worked
        if self.log > 1:
            print('LeucipPy(2) TEST IT WORKED #####')
            test1 = self.applyTransformations(transformations,round(coordsPole[0],4),round(coordsPole[1],4),round(coordsPole[2],4))
            print('LeucipPy(2) CoordsPole=',round(test1[0],4),round(test1[1],4),round(test1[2],4))
            test2 = self.applyTransformations(transformations,round(coordsCentre[0],4),round(coordsCentre[1],4),round(coordsCentre[2],4))
            print('LeucipPy(2) CoordsCentre=',round(test2[0],4),round(test2[1],4),round(test2[2],4))
            print('##### LeucipPy(2) TESTED #####')
        return transformations

    def applyTransformations(self,transforms,x,y,z):
        if self.log > 1:
            print('LeucipPy(2) Transform coords=',round(x,4),round(y,4),round(z,4))

        a_translate = transforms[0]
        x, y, z = x + a_translate[0], y + a_translate[1], z + a_translate[2]
        if self.log > 1:
            print('LeucipPy(2) A.translation=',a_translate,round(x,4),round(y,4),round(z,4))

        b_translate = transforms[1][0]
        x, y = self.rotate(x, y, b_translate)
        if self.log > 1:
            print('LeucipPy(2) B.rotation',b_translate,round(x,4),round(y,4),round(z,4))

        c_translate = transforms[2][0]
        x,z = self.rotate(x, z, c_translate)
        if self.log > 1:
            print('LeucipPy(2) C.rotation',c_translate,round(x,4),round(y,4),round(z,4))

        '''d_translate = transforms[3]
        x, y, z = x + d_translate[0], y + d_translate[1], z + d_translate[2]
        if self.log > 1:
            print('LeucipPy(2) D.translation',d_translate,round(x,4),round(y,4),round(z,4))'''

        return x,y,z

    def reverseTransformations(self,transforms,x,y,z):
        if self.log > 1:
            print('LeucipPy(2) Reverse coords=',round(x,4),round(y,4),round(z,4))

        '''d_translate = transforms[3]
        x, y, z = x - d_translate[0], y - d_translate[1], z - d_translate[2]
        if self.log > 1:
            print('LeucipPy(2) D.rev translation', d_translate,round(x,4),round(y,4),round(z,4))'''

        c_translate = transforms[2][0]
        x,z = self.rotate(x,z, 360-c_translate)
        if self.log > 1:
            print('LeucipPy(2) C.rev rotation',c_translate,round(x,4),round(y,4),round(z,4))

        b_translate = transforms[1][0]
        x, y = self.rotate(x, y, 360-b_translate)
        if self.log > 1:
            print('LeucipPy(2) B.rev rotation',b_translate,round(x,4),round(y,4),round(z,4))

        a_translate = transforms[0]
        x,y,z = x-a_translate[0],y-a_translate[1],z-a_translate[2]
        if self.log > 1:
            print('LeucipPy(2) A.rev translation=',a_translate,round(x,4),round(y,4),round(z,4))

        return x,y,z

    def rotate(self,x,y,degrees):
        '''
        Tis rotates the fiven coords around the orgin, leaving x-plane undisturbed
        :param x: the original position
        :param y: the original position
        :param degrees: how far to rotate in degrees
        :return: the new position  as a tuple x,y
        '''
        #print('This rotates an atom about the origin XY the given desgrees')
        radians = self.__degreesToRadians(degrees)
        angle_left = radians
        x_now, y_now = x,y
        while (angle_left > M_PI / 2):
            x_now,y_now = self.__rotateNinety(x_now, y_now)
            angle_left -= M_PI / 2
        x_now, y_now = self.__rotateQuadrant(x_now, y_now, angle_left)
        return x_now, y_now

    def makeRandomVersions(self,number):
        #print('This will make the given number of randomised geometry versions af the original as if new pdbs')
        # it will return a dataframe in a format that my LeuipPy library will process to do geometry
        pdbs = []
        resnos = []
        atoms = []
        aas = []
        xs = []
        ys = []
        zs = []
        for i in range(0,number):
            self.moveDataFrameToReadyLists()
            # first do all the shifting of the atoms to their new places
            for branch in self.branches_sets:
            #for branch, atom_list in self.branches_sets.items():
                if True:#branch == '1':
                    for j in range(1,len(self.branches_sets[branch])):
                        #atm_prev = self.branches_sets[branch][j-1]
                        atm_this = self.branches_sets[branch][j]
                        atm_prev = atm_this.last_atom
                        pole = atm_prev.atom + '-' + atm_this.atom
                        lx,ly,lz = atm_prev.coordX,atm_prev.coordY,atm_prev.coordZ
                        tx,ty,tz = atm_this.coordX,atm_this.coordY,atm_this.coordZ
                        trans = self.shiftToPole([lx,ly,lz],[tx,ty,tz])
                        degrees = self.branches_sets[branch][j].rotation.getRandomRotation()
                        #print(degrees)
                        #degrees = 360
                        if degrees != None:
                            if self.log > 0:
                                print('LeucipPy(1) Rotate degrees=', round(degrees, 4), 'around pole=', pole)
                            for k in range(j+1,len(self.branches_sets[branch])):
                                atm = self.branches_sets[branch][k]
                                nx,ny,nz = self.applyTransformations(trans,atm.coordX,atm.coordY,atm.coordZ)
                                if self.log > 1:
                                    print('LeucipPy(2) Rotate degrees=', round(degrees, 4), 'from=', round(nx, 4), round(ny, 4), round(nz, 4))
                                ny,nz = self.rotate(ny,nz,degrees)
                                if self.log > 1:
                                    print('LeucipPy(2) Rotate degrees=', round(degrees, 4), 'to=', round(nx, 4),round(ny, 4),round(nz, 4))
                                nx, ny, nz = self.reverseTransformations(trans,nx,ny,nz)
                                atm.coordX = nx
                                atm.coordY = ny
                                atm.coordZ = nz
                                self.branches_sets[branch][k] = atm




            # then all the adding, but we only add branch 1 as all atoms start on branch 1
            atoms_list = self.branches_sets['1']
            for atm in atoms_list:
                pdbs.append('pdb' + str(i))
                resnos.append(atm.resNo)
                atoms.append(atm.atom)
                aas.append(atm.aminoAcid)
                xs.append(atm.coordX)
                ys.append(atm.coordY)
                zs.append(atm.coordZ)

        dic_atoms = {}
        dic_atoms['pdb_code'] = pdbs
        dic_atoms['rid'] = resnos
        dic_atoms['aa'] = aas
        dic_atoms['atom'] = atoms
        dic_atoms['x'] = xs
        dic_atoms['y'] = ys
        dic_atoms['z'] = zs
        atoms_data = pd.DataFrame.from_dict(dic_atoms)
        return atoms_data

    def getGeometry(self,geos):
        print('This will get geometry like from my library but a bit different as it is not')
        print('Or it will return something that my library can process instead of biopython')
        print('Or it will cretae biopython objects that I can then send to my library')

    #######################################################################################################
    ### This code is copied from my C++ implentation in LeucipPus
    def __getRotationAngle(self, x, y):
        theta_radians = 0.0
        qStart = self.__getQuadrant(x, y)
        mag = math.sqrt(pow(x, 2) + pow(y, 2))
        newXY = 0,0
        axisXY = 0,0
        if (mag > 0.0001):
            newXY = x,y
            axisXY= mag,0
            theta_degrees = calc.getAngle(newXY[0],newXY[1],0,0,0,0,axisXY[0],axisXY[1],0)
            #theta_degrees = calc.getAngle(0, 0, 0,newXY[0], newXY[1], 0, axisXY[0], axisXY[1], 0)
            #theta_degrees = calc.getAngle(newXY[0], newXY[1], 0, axisXY[0], axisXY[1], 0,0, 0, 0)
            theta_radians = self.__degreesToRadians(theta_degrees)
            # We want to go claockwise to the  positive  x - axis and this is just the absolute difference, so:
            if (qStart == 4 or qStart == 3):
                radians = 2 * M_PI - theta_radians
            if (theta_radians < 0):
                theta_radians = 2 * M_PI + theta_radians
        return theta_radians

    def __getQuadrant(self,x,y):
        '''
        This gets toe quadrant on the given plane, note x,y could include z we are golding 1 plane fixed
        :param x:
        :param y:
        :return: the quadrant 1,2,3 or 4
        '''
        qStart = 1
        if (x == 0 and y > 0):
            qStart = 1
        elif (x == 0 and y < 0):
            qStart = 4
        elif (y == 0 and x > 0):
            qStart = 1
        elif (y == 0 and x < 0):
            qStart = 2
        elif (x < 0 and y < 0):
            qStart = 3
        elif (x < 0 and y > 0):
            qStart = 2
        elif (y < 0 and x > 0):
            qStart = 4
        return qStart

    def __rotateNinety(self,x,y):
        '''
        Rotates into another quadrant
        :param x:
        :param y:
        :return: x,y tuple
        '''
        q = self.__getQuadrant(x, y)
        q -= 1
        if (q == 0):
            q = 4
        next_x, next_y = abs(y), abs(x)
        if (q == 2):
            next_x *= -1
        elif (q == 3):
            next_x *= -1
            next_y *= -1
        elif (q == 4):
            next_y *= -1
        return next_x,next_y

    def __rotateQuadrant(self,x,y,radians):
        if (abs(radians) > 0.001):
            # This assumes  an  angle that is positive and less or = than   90  that   may  turn  only  into   the  next quadrant.
            # Choose quadrants
            qStart = self.__getQuadrant(x, y)
            mag = math.sqrt(math.pow(x, 2) + math.pow(y, 2)) # the length of the vector
            if (mag > 0.0001):
                sinA = abs(y) / mag
                angleA = math.asin(sinA) # this is the    angle    made   with the x-axis from the original vector
                angleB = angleA - radians # this is the angle made with the x-axis with the rotated vector
                qEnd = qStart
                if (qStart == 1):
                    if (radians > angleA):
                        angleB = radians - angleA
                        qEnd = 4
                elif (qStart == 2):
                    angleB = radians + angleA
                    if (angleA + radians > M_PI / 2):
                        angleB = M_PI - (angleA + radians)
                        qEnd = 1
                elif (qStart == 3):
                    if (radians > angleA):
                        angleB = radians - angleA
                        qEnd = 2
                else: # must be q4
                    angleB = radians + angleA
                    if (angleA + radians > M_PI / 2):
                        angleB = M_PI - (radians + angleA)
                        qEnd = 3
                x2 = math.cos(angleB) * mag
                y2 = math.sin(angleB) * mag
                if (qEnd == 2 or qEnd == 3):
                    x2 *= -1
                if (qEnd == 3 or qEnd == 4):
                    y2 *= -1
                return x2,y2
            else:
                return x,y
        else:
            return x,y

    def __degreesToRadians(self,degrees):
        return (degrees * M_PI)/180
    def __radiansToDegrees(self,radians):
        return (radians * 180)/M_PI







