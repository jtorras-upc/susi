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

import re
import os           # Allow use of function path.exists
######################################################################################################################
class parserSusi():
    __version__ = '2.0.0'
######################################################################################################################
    def __init__(self):
        # Assigning default values to INPUT variables
        self.box                              = []
        self.periodical_coordinates           = 'N'
        self.minimum_distance                 = [0.8,'A']
        self.retries                          = 5
        self.maximum_deviation                = 5.0
        self.t_errors                         = []
        self.t_molecules                      = []
        self.systemName                       = ''
        self.pathPrepi                        = ''
        self.t_crosslinks                     = []
        self.t_starting_closing_atoms         = []
        self.t_sub_starting_closing_atoms     = []
        self.molecules_number                 = 0
        self.crosslinks_number                = 0
        self.bond_length_start                = [-1,'A']
        self.bond_length_close                = [-1,'A']
        self.distorion_angle                  = -1
        self.t_cyclic                         = []
        
        # Dictionary with all regular expressions to be used in the Input file parsing
        self.rx_dict = {
             'system_name': re.compile(r'system_name(.*,)(?P<system_name>.*)\n'),
             'box': re.compile(r'box(.*,)Lx(?P<Lx>.*)(?P<Ux>.)(.*,.*)Ly(?P<Ly>.*)(?P<Uy>.)(.*,.*)Lz(?P<Lz>.*)(?P<Uz>.)\n'),
             'periodical_coordinates': re.compile(r'periodical_coordinates(.*,)(?P<Pc>.*)\n'),
             'minimum_distance': re.compile(r'minimum_distance(.*,)(?P<Md>.*)(?P<U>.)\n'),
             'retries': re.compile(r'retries(.*,)(?P<retries>.*)\n'),
             'maximum_deviation': re.compile(r'maximum_deviation(.*,)(?P<MaxD>.*)\n'),
             'molecules': re.compile(r'molecules(.*,)(?P<MolN>.*)\n'),
             'path_prepi': re.compile(r'path_prepi(.*,)(?P<PathP>.*)\n'),
             'crosslinks': re.compile(r'crosslinks(.*,)(?P<CLN>.*)\n'),
             'bond_length_CL_start': re.compile(r'bond_length_CL_start(.*,)(?P<BLCLS>.*)(?P<Ub>.)\n'),
             'bond_length_CL_close': re.compile(r'bond_length_CL_close(.*,)(?P<BLCLC>.*)(?P<Ub>.)\n'),
             'distortion_angle': re.compile(r'distortion_angle(.*,)(?P<DAng>.*)\n'),
        }

######################################################################################################################
######################################################################################################################
    #@staticmethod
    def getAndCheckFileValues(self,t_lines):
        """Reads the information in the inputfile and check values.

        Input-----------------------
            t_lines (str): Lines from the input file.
        Output----------------------  
            box (list): List with 3 box dimensions (X,Y,Z) and units.
            periodical_coordinates (str): Use boundary condition / N: Don't use boundary condition.
            minimum_distance (list): Size and unit of minimum distance to allow particle collisions.
            retries (int): Maximum number of retries to find a residue angle that don't collide with other atoms
            maximum_deviation (float): Value for maximum deviation.
            systemName (str): System name.
            pathPrepi (str): Path of PREPI files.
            t_molecules (list): List of names of molecules obtained from the input file.
            t_crosslinks (list): List of names of crosslink molecules obtained from the input file.
            t_starting_closing_atoms (list): List with the starting and closing atoms names.
            t_sub_starting_closing_atoms (list): List with the atoms that can be subsituted.
            distorion_angle (float): value of the distortion angle.
            bond_length_start (list): bond length and unit for the starting atom.
            bond_length_close (list): bond length and unit for the closing atom.
            t_errors (list): If some errors occured, this list contain all of them
        """

        # Boolean variables of mandatory input values
        isSystemName      = False
        isBox             = False
        isPath            = False
        isMolecules       = False
        isCrosslinks      = False
        isAngle           = False
        isBondLenghtStart = False
        isBondLenghtClose = False
        # Auxiliary counters
        auxMol            = 0
        auxCL             = 0

        # Parsing INPUT file
        ####################
        for line in t_lines:
            key, match = self.parse_line(line)

            if key == "system_name":
                self.systemName = match.group('system_name').strip()
                isSystemName = True
            elif key == "box":
                self.box.append(float(match.group('Lx').strip()))  # X length
                self.box.append(      match.group('Ux').strip())   # X Units
                self.box.append(float(match.group('Ly').strip()))  # Y length
                self.box.append(      match.group('Uy').strip())   # Y Units
                self.box.append(float(match.group('Lz').strip()))  # Z length
                self.box.append(      match.group('Uz').strip())   # Z Units
                isBox = True
            elif key == "periodical_coordinates":
                self.periodical_coordinates = match.group('Pc').strip()
            elif key == "minimum_distance":
                self.minimum_distance[0] = float(match.group('Md').strip())                
                self.minimum_distance[1] = match.group('U').strip()
            elif key == "retries":
                self.retries = int(match.group('retries').strip())
            elif key == "maximum_deviation":
                self.maximum_deviation = float(match.group('MaxD').strip())
            elif key == "molecules":
                self.molecules_number = int(match.group('MolN').strip())
                isMolecules = True
            elif key == "path_prepi":
                self.pathPrepi = match.group('PathP').strip()
                isPath = True
            elif key == "crosslinks":
                self.crosslinks_number = int(match.group('CLN').strip())
            elif key == "bond_length_CL_start":
                self.bond_length_start[0] = float(match.group('BLCLS').strip())
                self.bond_length_start[1] = match.group('Ub').strip()
                isBondLenghtStart = True 
            elif key == "bond_length_CL_close":
                self.bond_length_close[0] = float(match.group('BLCLC').strip())
                self.bond_length_close[1] = match.group('Ub').strip()
                isBondLenghtClose = True            
            elif key == "distortion_angle":
                self.distorion_angle = float(match.group('DAng').strip())
                isAngle = True
                isCrosslinks = True

            # Get all molecules: PREPI files
            elif isMolecules and not isCrosslinks :
                auxMol += 1
               # Check if more than one residue is reported. Ej.: 5 RES1
                if len(line.strip()) > 0 :   # Ignore blank lines
                  # Get all values of this line separated by ',' an store in list t_aux        
                    if line.startswith("[c]"): #check if there is a flag for cyclic polymers
                        t_aux = line[3:].split(",")
                        t_aux = [i.strip() for i in t_aux] # Remove string '\n' from end of each line
                        t_aux = parserSusi.checkRepeatedResidues(t_aux)
                        self.t_molecules.append(t_aux)
                        self.t_cyclic.append(True)
                    else: 
                        t_aux = line.split(",") 
                        t_aux = [i.strip() for i in t_aux] # Remove string '\n' from end of each line
                        t_aux = parserSusi.checkRepeatedResidues(t_aux)
                        self.t_molecules.append(t_aux)
                        self.t_cyclic.append(False)
                  # Closing molecule parsing
                    if auxMol >= self.molecules_number:  isMolecules = False

            # Get all crosslinks: PREPI files
            elif isCrosslinks and isBondLenghtStart and isBondLenghtClose and isAngle:
                auxCL += 1
                #Extract starting and closing atoms names
                try:    
                    #Starting atom
                    start_index = line.index('[') + 1
                    end_index = line.index(':')
                    start_atom, sub_start_atom = line[start_index:end_index].split(',')

                    #Atom to remove from the starting atom residue                                        
                    #Closure atom
                    start_index = end_index + 1
                    end_index = line.index(']')
                    closing_atom, sub_closing_atom = line[start_index:end_index].split(',')

                    self.t_starting_closing_atoms. append([start_atom.strip(),closing_atom.strip()])
                    self.t_sub_starting_closing_atoms. append([sub_start_atom.strip(),sub_closing_atom.strip()])

                    # Extract the list of residues
                    start_index = line.index(']') + 1
                    end_index = len(line)
                    line_aux = line[start_index:end_index]
                # Check if more than one residue is reported. Ej.: 5 RES1
                    if len(line_aux.strip()) > 0 :   # Ignore blank lines
                    # Get all values of this line separated by ',' an store in list t_aux
                        t_aux = line_aux.split(",") 
                        t_aux = [i.strip() for i in t_aux] # Remove string '\n' from end of each line
                        t_aux = parserSusi.checkRepeatedResidues(t_aux)
                        self.t_crosslinks.append(t_aux)
                    # Closing crosslink parsing
                        if auxCL >= self.crosslinks_number:  isCrosslinks = False
                except ValueError:    # Handle the case where [] is not found
                    self.t_starting_closing_atoms = [[None]]
                    self.t_sub_starting_closing_atoms = [[None]]

         # Checking for errors and/or missing values
         ###########################################
        if not isSystemName: self.t_errors.append('003') #system_name not found in file
        if not isBox       : self.t_errors.append('004') #box dimension not found
        if not isPath      : self.t_errors.append('010') #PATH of prepi files not found in file

        if len(self.periodical_coordinates)== 0: 
            self.t_errors.append('005')  #periodical_coordinates not found
        if len(self.minimum_distance)      == 0:
            self.t_errors.append('006')  #minimum_distance not found
        if self.retries                    <= 0:
            self.t_errors.append('007')  #retries not found in file
        if self.maximum_deviation          <= 0.0:
            self.t_errors.append('008')  #maximum_deviation not found in file
        if self.molecules_number                == 0: 
            self.t_errors.append('009')  #molecules_number not found in file

        if len(self.box) < 6:               self.t_errors.append('022')        #Error in box

         # Check units type. At this moment only Amstrong's are supported
        if self.box[1] != 'A':     self.t_errors.append('011')   #Error in X dimension unit
        if self.box[3] != 'A':     self.t_errors.append('012')   #Error in Y dimension unit
        if self.box[5] != 'A':     self.t_errors.append('014')   #Error in Z dimension unit
        if self.minimum_distance[1] !='A':
            self.t_errors.append('015')   #Error in minimum distance unit

         #System name not found in file
        if self.systemName == '':  self.t_errors.append('016')

         #Error in value of periodical coordinates. Possible values are Y or N
        if self.periodical_coordinates != 'Y' and self.periodical_coordinates != 'N':
            self.t_errors.append('017')

         #Error in minimum distance   
        if len(self.minimum_distance) < 2: self.t_errors.append('018')

         #Error in minimum distance units  
        if self.minimum_distance[1] != 'A':
            self.t_errors.append('035')    

         #Molecules number not found in file
        if self.molecules_number == 0:  self.t_errors.append('019') 

         #PATH of prepi files not found in file
        if self.pathPrepi   =="":
            self.t_errors.append('020')
        elif os.path.isdir(self.pathPrepi) == "False": 
            #self.t_errors.append("Could not found PATH:"+pathPrepi)
            self.t_errors.append("Could not found PATH")

         #Molecules number in file different than molecules number found
        if int(self.molecules_number) != len(self.t_molecules):
            self.t_errors.append('021')   

        #Wrong value for ditortion angle
        if isAngle and self.distorion_angle <= 0: 
            self.t_errors.append('029')   
        
        if isCrosslinks and isAngle == False: 
            self.t_errors.append('029')   
        
        #Wrong value for bond length
        if isBondLenghtStart and self.bond_length_start[0] <= 0: 
            self.t_errors.append('034')   
        
        if isCrosslinks and isBondLenghtStart == False: 
            self.t_errors.append('034')   
        
        #Error in the bond length units
        if isBondLenghtStart and self.bond_length_start[1] != 'A':
            self.t_errors.append('035')    

        #Error in minimum distance greater than Crosslink bond length
        if isBondLenghtStart and self.bond_length_start[0] < self.minimum_distance[0]:
            self.t_errors.append('036')    
                
        #####
        #Wrong value for bond length
        if isBondLenghtClose and self.bond_length_close[0] <= 0: 
            self.t_errors.append('041')   
        
        if isCrosslinks and isBondLenghtClose == False: 
            self.t_errors.append('041')   
        
        #Error in the bond length units
        if isBondLenghtClose and self.bond_length_close[1] != 'A':
            self.t_errors.append('035')    

        #Error in minimum distance greater than Crosslink bond length
        if isBondLenghtClose and self.bond_length_close[0] < self.minimum_distance[0]:
            self.t_errors.append('042')    

        #Warning minimum distance greater than C-C single bond length
        if self.minimum_distance[0] > 1.54:
            print('Warning: Minimum distance is greater than that of a typical C-C single bond'+'\n')
    
        #Warning crosslink bond length is greater than a C-C single bond length
        if isBondLenghtStart and self.bond_length_start[0] > 1.54:
            print('Warning: Starting Crosslink bond length is greater than that of a typical C-C single bond'+'\n')

        #Warning crosslink bond length is greater than a C-C single bond length
        if isBondLenghtClose and self.bond_length_close[0] > 1.54:
            print('Warning: Closing Crosslink bond length is greater than that of a typical C-C single bond'+'\n')

        #Starting and closing atoms not found
        for pair_atoms in self.t_starting_closing_atoms:
            if len(pair_atoms) != 2 or pair_atoms[0] == None or pair_atoms[1] == None or self.t_starting_closing_atoms == [] :
                self.t_errors.append('027')   
            
        #Subtitution Starting and closing atoms not found
        for pair_atoms in self.t_sub_starting_closing_atoms:
            if len(pair_atoms) != 2 or pair_atoms[0] == None or pair_atoms[1] == None or self.t_sub_starting_closing_atoms == [] :
                self.t_errors.append('038')   

         #Crosslinks number in file different than crosslinks number found
        if self.crosslinks_number != len(self.t_crosslinks):
            self.t_errors.append('031')   
        
        return self.box, self.periodical_coordinates, self.minimum_distance, self.retries,  \
                  self.maximum_deviation, self.systemName, self.pathPrepi,                      \
                  self.t_molecules, self.t_crosslinks, self.t_starting_closing_atoms, self.t_sub_starting_closing_atoms, self.distorion_angle, self.bond_length_start, self.bond_length_close, self.t_cyclic, self.t_errors

######################################################################################################################
######################################################################################################################
    @staticmethod
    def checkRepeatedResidues(t_residues):
        """Reads all residue list from a molecule and checks if some of this residue values indicates than this residue has to be repeated.
        Input-----------------------
        t_residues (list): Residue list from a molecule name.
        Output----------------------
        t_auxResidues (list): Returns residue list with repeated values in position with number and residue name.
        """
        t_auxResidues = []
        for v_lines_number in range(0,len(t_residues)):
            t_residue = t_residues[v_lines_number].split(" ")
            if len(t_residue) == 2:
                  # More than one value in residue separated by ','. Example 3 ALA
                v_numEqualResidues = int(t_residue[0])  #Number of same residues
                v_counter = 0
                while v_numEqualResidues > v_counter:   #Add number times to t_molecules
                    t_auxResidues.append(t_residue[1])
                    v_counter +=1
            else:
                  # Only one value in residue. Example ALA
                t_auxResidues.append(t_residues[v_lines_number])
        return t_auxResidues

######################################################################################################################
######################################################################################################################
    def parse_line(self,line):
        """Reads a line as a string and look for a match on the dictionary of regular expressions to parse an input file.

        Input-----------------------
            line (str):  Line to check for matches of regular expressions
        Output----------------------
            t_aux (list): Dictionary key for the regular expression found in the line
            match (list): Regular expression object found in the linke
        """
        for key, rx in self.rx_dict.items():
            match = rx.search(line)
            if match:
             #print "FOUNDED ==>>",key
               return key,match
 
        return None, None
######################################################################################################################
