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

import sys     # To allow use of function sys.exit()
import os.path # To get file extensions with os.path.splitext

class filePREPI:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self):
        self.name         = ''          #Residue name
        self.t_file_lines = []
        self.t_name       = []          #Atom name
        self.t_con        = []
        self.t_par        = []
        self.t_mLines     = []          #List with lines with M values
        self.t_topologic  = []          #List with lines with topologycal type
        self.t_numbering  = []          #List with lines with the numbering of the atoms 
######################################################################################################################
######################################################################################################################
    def readFile(self,fileName):
        """Check if file exists, if not exists returns a error message

        Input-----------------------
            fileName (str): Name of current file.
        """
        try:
            v_file     = open(fileName, "r")
            self.name  = os.path.splitext(os.path.basename(fileName))[0]
            file_lines = []
    # Store all PREPI lines in file_lines list
            for v_line in v_file:
                file_lines.append(v_line)
            v_file.close()
            v_read_file = ""
            for v_num_lines in range(0,len(file_lines)):
                if file_lines[v_num_lines][6:10] == "DUMM" and v_read_file == "":
                    v_read_file = "X"
                    self.t_file_lines.append(file_lines[v_num_lines])
                else:
                    if v_read_file == "X":
                        if file_lines[v_num_lines][0:4] == "\n":
                            break                # No more lines from file have to be readed
                        else:
    #                       if file_lines[v_num_lines][6:10] <> "DUMM":
                            self.t_file_lines.append(file_lines[v_num_lines])
            self.__trataFile__()
        except IOError:
# File not exists, a error message would be displayer and exit of current run
            print("Could not read file:", fileName)
            sys.exit()
######################################################################################################################
######################################################################################################################
    def __trataFile__(self):
        """Fill values for t_name, t_con, t_par, that would be used in int2car_M class.
        """

        v_num_lines = 0
        while v_num_lines < len(self.t_file_lines):
            v_line = linePREPI(self.t_file_lines[v_num_lines])
            self.t_name.append(v_line.igraph)
            self.t_con.append(filePREPI.getCar(v_line.na, v_line.nb, v_line.nc))
            self.t_par.append(filePREPI.getPar(v_line.r, v_line.theta, v_line.phi))
            self.t_topologic.append(v_line.itree) #######################Shows the type of atom (M, B, E, etc)
            self.t_numbering.append(v_line.number) #######################Shows numbering of the atoms in the PREPI file
            if v_line.itree == 'M':
                self.t_mLines.append(v_num_lines)
            v_num_lines +=1
######################################################################################################################
######################################################################################################################
    @staticmethod
    def getCar(na, nb, nc):
        """Create a list with values of the connectiviy of atoms according to PREPI format.

        Input-----------------------
            na (str): atom index for bond reference.
            nb (str): atom index for angle reference.
            nc (str): atom index for dihedral reference.
        Output----------------------
            t_aux (list): list with values of the connectiviy of atoms according to PREPI format.
        """
        t_aux = []
        t_aux.append(na)
        t_aux.append(nb)
        t_aux.append(nc)
        return t_aux
######################################################################################################################
######################################################################################################################
    @staticmethod
    def getPar(r, theta, phi):
        """Create a list with values of the internal coordinates according to PREPI format.

        Input-----------------------
            r (str): bond length.
            theta (str): planar angle.
            phi (str): dihedral angle.
        Output----------------------
            t_aux (list): List with values of r, theta and phi as float.
        """

        t_aux = []
        t_aux.append(float(r))
        t_aux.append(float(theta))
        t_aux.append(float(phi))
        return t_aux
######################################################################################################################
######################################################################################################################
    def checkFile(self):
        """Print current file.
        """
        print("Residue name:  " + str(self.name))
        print("Atoms name:    " + str(self.t_name))
        print("Con:           " + str(self.t_con))
        print("Par:           " + str(self.t_par))
        print("M lines:       " + str(self.t_mLines))
######################################################################################################################
######################################################################################################################
    @staticmethod
    def storeFile():
        """Empty list to be created later if it's necessary
        """
        print("PREPI file can not be generated")
######################################################################################################################
######################################################################################################################
class linePREPI:
    def __init__(self, line):
        """Class with all the information extracted from a PREPI file.

        Input-----------------------
        line (str): line from a PREPI file.

        """
        self.line                    = line
        self.number                  = int(line[0:4])      # I: Current number of the atom in the tree
        self.igraph                  = line[6:10].strip()  # Maximum 4 characters
        self.isymbl                  = line[12:14]
        self.itree                   = line[18:19]         # Topological type (M, S, B, E, or 3)
        self.na                      = int(line[21:24])    # Atom number to which atom I is connected
        self.nb                      = int(line[25:28])
        self.nc                      = int(line[29:32])
        self.r                       = line[37:42].strip() # coordinates
        self.theta                   = line[45:52].strip() # coordinates
        self.phi                     = line[54:62].strip() # coordinates
        self.chrg                    = line[63:72]         # Partial atomic charge on atom I
######################################################################################################################
