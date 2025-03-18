#########################################################################
#
#                          ---  SuSi   ---                          
#
#                       Structure Simulation                        
#                                                                   
#           Copyright (C) 2025                                      
#                                                                   
#           David Naranjo        david.alejandro.naranjo@upc.edu    
#           Carlos Alemán        carlos.aleman@upc.edu              
#           Juan Torras          joan.torras@upc.edu                
#                                                                   
#########################################################################
#     This file is part of SuSi.
#
# SuSi is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your 
# option) any later version.
#
# SuSi is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
# for more details.
#
# You should have received a copy of the GNU General Public License 
# along with SuSi. If not, see <https://www.gnu.org/licenses/>. 
#
# Neither the names of the EEBE, the Universitat Politecnica de 
# Catalunya, nor the names of any of the copyright holders of SuSi may 
# be used to endorse or promote any products derived from this Software 
# without specific, prior, written permission. 
#########################################################################

import math
import random
import numpy as np
import sys

class int2car:
    global int2carx, int2cary, int2carz, new_car
    """this program tranformates from internal to cartesian coordinates. This is done with the NERF(Natural Extension Reference Frame) method
    NERF has 2 steps, first the position of the D(new) atom is given with the spherical coordinate representation with C as a reference
    then the ABC reference frame is introduced an D recalculated to fit this frame
    """
########################################################
########################################################
    def int2carx(dist,alpha,dihed):
        """Returns the value of the X-axis coordinate from internal coordinates.

        Input-----------------------
            dist (float): bond length.
            alpha (float): planar angle.
            dihed (float): dihedral angle.
        Output----------------------
            x (float): X-axis coordinate.
        """
        alpha = math.radians(alpha)
        #dihed = math.radians(dihed) 
        x =  dist*math.cos(math.pi-alpha) 	
        return x
########################################################
########################################################
    def int2cary(dist,alpha,dihed):
        """Returns the value of the Y-axis coordinate from internal coordinates.

        Input-----------------------
            dist (float): bond length.
            alpha (float): planar angle.
            dihed (float): dihedral angle.
        Output----------------------
            y (float): Y-axis coordinate.
        """
        alpha = math.radians(alpha)
        dihed = math.radians(dihed) 
        y =  dist*math.cos(dihed)*math.sin(math.pi-alpha)
        return y
########################################################
########################################################
    def int2carz(dist,alpha,dihed):
        """Returns the value of the Z-axis coordinate from internal coordinates.

        Input-----------------------
            dist (float): bond length.
            alpha (float): planar angle.
            dihed (float): dihedral angle.
        Output----------------------
            z (float): Z-axis coordinate.
        """
        alpha = math.radians(alpha)
        dihed = math.radians(dihed) 
        z =  dist*math.sin(dihed)*math.sin(math.pi-alpha)
        return z
########################################################
########################################################    
    def new_car(dist,alpha,dihed,pC,pB,pA):
        """Gets a new cartesian coordinate from the new reference frame as a funtion of the 3 last atoms(A,B,C).

        Input-----------------------
            dist (float): bond length.
            alpha (float): planar angle.
            dihed (float): dihedral angle.
            pC (list): coordinates of the C atom.
            pB (list): coordinates of the B atom.
            pA (list): coordinates of the A atom.
        Output-----------------------
            newcoord (numpy.ndarray): new cartesian coordinates from internal.
        """
        C = pC  #coordinates of the C atom with relative distance dist
        coordinatex = int2carx(dist,alpha,dihed)
        coordinatey = int2cary(dist,alpha,dihed)
        coordinatez = int2carz(dist,alpha,dihed)
        D2 = np.array([coordinatex, coordinatey, coordinatez]) #array with the displacement with C as reference
        AB = (pB - pA) #parameters neeeded for the reference frame transformation
        BC = (pC - pB) #pB are the coordinates for the atom that forms the angle current,pC,pB. dA is dihedral analogously
        bc = BC/np.linalg.norm(BC)
        n = np.cross(AB,bc)/np.linalg.norm(np.cross(AB,bc))
        M = np.array([bc,np.cross(n,bc),n])   #array of arrays(matrix) with the conversion to the ABC reference frame  
        newcoord = C + np.dot(D2,M)       #tranformation from internal to cartesian
        
        return newcoord
########################################################
    @staticmethod
    def print_xyz(title, names, coord):
        """Prints the cartesian coordinates.

        Input-----------------------
            title (str): title of the message.
            names (list): list of atoms.
            coord (list): list of coordinates.
        """

        n = len(names)
        print(n)
        print(title)
        for i in range(n):
            print(names[i][0],coord[i][0],coord[i][1],coord[i][2])

########################################################
    @staticmethod
    def int2car(con,par,xyz):
        """Computes the cartesian coordinates.

        Input-----------------------
        con (list): list of lists of integers representing the connectivity of atoms according to PREP format.
        par (list): list of lists of floats representing the parameters of internal coordinates of atoms according to PREP format
        xyz (list): list of lists of floats representig the cartesian coordinates.

        Output----------------------
        xyz (list): list of lists of floats representig the cartesian coordinates.
        """
        #if xyz list is given empty then the first 3 atoms are putted as follows
        #   pdb.set_trace()
        if len(xyz) == 1 or len(xyz) == 2:
        # xyz list has to be zero or 3 or more elements 
            print("Bad numbers of z-matrix elements",len(xyz)," (should be 0,3 or more):")
            sys.exit()
        elif len(xyz)==0:      
        #1st atom is the origin
            pC = np.array([0.0,0.0,0.0])  
            xyz.append(pC)
        #2nd atom will bi situated randomly in the sapece at a given distance par[1][0]
            random.seed     # Create a completely new random serie
            rphi = random.uniform(0, 2*math.pi)
            rpsi = random.uniform(0, 2*math.pi) 
            xyz.append(np.array([par[1][0]*math.sin(rphi)*math.cos(rpsi),
            par[1][0]*math.sin(rphi)*math.sin(rpsi),
            par[1][0]*math.cos(rphi)]))
        #3rd atom is in the plane of previous 2 particles at a known distance of the 
        #2nd atom and forming a known angle with the first atom
            if (con[2][0])== 2:
                ang = math.pi - math.radians(par[2][1])
            else:
                ang = math.radians(par[2][1])
            xyz.append(np.array([xyz[con[2][0]-1][0] + par[2][0]*math.cos(ang),
            xyz[con[2][0]-1][1] + par[2][0]*math.sin(ang),
            xyz[con[2][0]-1][2]]))

        i = 3
        while i<len(con):
            coord = new_car(par[i][0],par[i][1],par[i][2],xyz[con[i][0]-1],xyz[con[i][1]-1],xyz[con[i][2]-1])
            xyz.append(coord)
            i = i + 1
        return xyz
########################################################
