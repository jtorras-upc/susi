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

import os.path
import numpy as np  # Allow functions to calculate remainder between 2 numbers
import sys          # To allow use of function sys.exit()
import string
from .atom     import atom
from .molecule import molecule
from .residue  import residue
from .system   import system
from .messages import messages
class filePDB:
    __version__ = '2.0.0'
    t_file_lines = []
######################################################################################################################
######################################################################################################################
    def readFile(self,fileName,systemName):

        """Reads a PDB file.

        Input-----------------------
            fileName (str): Name of current file.
            systemName (str): Name of current system.
        Output----------------------
            v_system (system): System object.
        """
        state=[0,0,0]
        try:
            v_file       = open(fileName, "r")
            for v_line in v_file:                               # Read all the file, line by line
                if v_line[0:6] in ["REMARK"]:
                    #remark_line=v_line[11:-1]#end instead of -1 eas at the beginning
                    remark_line=v_line[7:-1].strip()
                    temp=[int(x) for x in remark_line.split(' ') ]
                    state[0]=temp[0]
                    state[1]=temp[1]
                if v_line[0:6] in ["ATOM  "]:                    # Select only ATOM ines
                    self.t_file_lines.append(v_line)             # Add lines of type PDB_line to process later
            v_file.close()
        except IOError:
            print("Could not read file:", fileName)
            sys.exit()
###################################
        if len(self.t_file_lines) == 0:
            print("No PDB lines found in file:", fileName)
            sys.exit()                                      # If no ATOM or HETATM line found, leave function
###################################
        v_system = system(systemName,0,0,20,[45,"A",45,"A",45,"A"],0,0)                      # Init system var, later will be return value of class
        v_system.state=state
# Init variables for creation of new residue for this molecule
        v_molecule                = molecule(fileName, self.t_file_lines[0][21:22])
        v_residue_name            = self.t_file_lines[0][17:20]             # Get first value of residue name
        t_atom                    = []                                      # Atom list for create residue
        v_residue_sequence_number = self.t_file_lines[0][22:26]             # Get first value of residue sequence number
        v_num_lines               = 0                                       # Variable num linia PDB

        while v_num_lines < len(self.t_file_lines ):                        # Read from first to last position
            if v_molecule.ID != self.t_file_lines[v_num_lines][21:22]:        # Add new molecule
                v_molecule.addResidue(residue(t_atom, v_residue_name, v_residue_sequence_number,v_molecule.ID))
                v_system.addMolecule(v_molecule)
                v_molecule                = None
                v_molecule                = molecule(fileName, self.t_file_lines[v_num_lines][21:22])
                v_residue_name            = self.t_file_lines[v_num_lines][17:20]
                t_atom                    = []                       # Delete line contents
                v_residue_sequence_number = self.t_file_lines[v_num_lines][22:26]
    #           t_atom.append(atom(t_file_lines[v_num_lines].line))  # Add new atom to future residue
                t_atom.append(atom( self.t_file_lines[v_num_lines][0:6],              #ATOM o HETATM
                                    int(self.t_file_lines[v_num_lines][6:11]),        #serial number
                                    self.t_file_lines[v_num_lines][12:16].strip(),    #name
                                    self.t_file_lines[v_num_lines][21:22],            #chain ID
                                    self.t_file_lines[v_num_lines][30:38],            #X
                                    self.t_file_lines[v_num_lines][38:46],            #Y
                                    self.t_file_lines[v_num_lines][46:54],            #Z
                                    self.t_file_lines[v_num_lines][54:60],            #occupancy
                                    self.t_file_lines[v_num_lines][60:66],            #temperature
                                    self.t_file_lines[v_num_lines][76:78],            #element_symbol
                                    self.t_file_lines[v_num_lines][78:80],            #charge
                                    self.t_file_lines[v_num_lines][17:20].strip(),    #residue_name
                                    int(self.t_file_lines[v_num_lines][22:26]),
                                    "E",1,[1,2,3]))      #residue_number
            else:
# Check if this residue is the same as previous record
                if v_residue_sequence_number == self.t_file_lines[v_num_lines][22:26]:
                    v_atom = atom( self.t_file_lines[v_num_lines][0:6],              #ATOM o HETATM
                                int(self.t_file_lines[v_num_lines][6:11]),
                                self.t_file_lines[v_num_lines][12:16].strip(),    #Name
                                self.t_file_lines[v_num_lines][21:22],            #chain ID
                                self.t_file_lines[v_num_lines][30:38],            #X
                                self.t_file_lines[v_num_lines][38:46],            #Y
                                self.t_file_lines[v_num_lines][46:54],            #Z
                                self.t_file_lines[v_num_lines][54:60],            #occupancy
                                self.t_file_lines[v_num_lines][60:66],            #temperature
                                self.t_file_lines[v_num_lines][76:78],            #element_symbol
                                self.t_file_lines[v_num_lines][78:80],            #charge
                                self.t_file_lines[v_num_lines][17:20].strip(),    #residue_name
                                int(self.t_file_lines[v_num_lines][22:26]),       #residue_number
                                "E",1,[1,2,3])                      
                    t_atom.append(v_atom) # Add new atom to the next residue that will be created
                else:                    # Residue is different than previous record
# Add new residue with atoms value
                    v_molecule.addResidue(residue(t_atom, v_residue_name, v_residue_sequence_number,v_molecule.ID))
# Init variables
                    v_residue_name            = self.t_file_lines[v_num_lines][17:20]
                    v_residue_sequence_number = self.t_file_lines[v_num_lines][22:26]
                    t_atom                    = []                                         # Delete list value
        #              v_atom                    = atom(t_file_lines[v_num_lines].line)
                    v_atom = atom( self.t_file_lines[v_num_lines][0:6],              #ATOM o HETATM
                                    int(self.t_file_lines[v_num_lines][6:11]),
                                    self.t_file_lines[v_num_lines][12:16].strip(),    #Name
                                    self.t_file_lines[v_num_lines][21:22],            #chain ID
                                    self.t_file_lines[v_num_lines][30:38],            #X
                                    self.t_file_lines[v_num_lines][38:46],            #Y
                                    self.t_file_lines[v_num_lines][46:54],            #Z
                                    self.t_file_lines[v_num_lines][54:60],            #occupancy
                                    self.t_file_lines[v_num_lines][60:66],            #temperature
                                    self.t_file_lines[v_num_lines][76:78],            #element_symbol
                                    self.t_file_lines[v_num_lines][78:80],            #charge
                                    self.t_file_lines[v_num_lines][17:20].strip(),    #residue_name
                                    int(self.t_file_lines[v_num_lines][22:26]),       #residue_number
                                    "E",1,[1,2,3])
                    t_atom.append(v_atom)     # Add new atom to the next residue that will be created
            v_num_lines += 1

        if (v_molecule.ID == self.t_file_lines[v_num_lines-1][21:22]):
            v_molecule.addResidue(residue(t_atom, v_residue_name, v_residue_sequence_number,v_molecule.ID))
            v_system.addMolecule(v_molecule)
        return v_system
######################################################################################################################
######################################################################################################################
    @staticmethod
    def storeFile(fileName,system,language):
        """Stores current file.

        Input-----------------------
        fileName (str): Name of current file.
        system (system): System object.
        """

        v_serial_number = 0
# Checks if exists a file with same name
        print("NUMMOLS to print",len(system.t_molecules))
        if os.path.isfile(fileName) == True:
            print('File', fileName, 'exists')
        else:    #File does not exists
            v_file = open(fileName, "w")

            #Write cell info

            Lx = f"{float(system.box[0]):9.3f}"
            Ly = f"{float(system.box[2]):9.3f}"
            Lz = f"{float(system.box[4]):9.3f}"

            alpha = beta = gamma = f"{float(90):7.3f}"

            Line_cell = f"CRYST1{Lx}{Ly}{Lz}{alpha}{beta}{gamma} P 1        1\n"
            v_file.write(Line_cell)

            if system.state[2]==1:
                v_file.write('REMARK   ' + str(system.state[0]) + ' ' + str(system.state[1]) +  str(system.state[2]) + ' \n') 
            for v_num_molecule in range(0,len(system.t_molecules)):
                for v_num_residue in range(0,len(system.t_molecules[v_num_molecule].t_residues)):
                    for v_num_atom in range(0,len(system.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom)):
    # no store in file DUMM lines
                        if system.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom[v_num_atom].name != 'DUMM':
                            v_serial_number +=1
    #                       line = system.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom[v_num_atom].addLine(v_serial_number)
                            line = filePDB.addLine(system.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom[v_num_atom], v_serial_number,language)
                            v_file.write(line)
    # Add end line
                #v_serial_number +=1

                line = filePDB.addEndLine2()

                #if system.state[0] > v_num_molecule:
                #   print system.state[0], v_num_molecule
                #   line = filePDB.addEndLine(system.t_molecules[v_num_molecule].t_residues[-1].t_atom[-1], v_serial_number)
                #else:
                #   print system.state[0], v_num_molecule
                #   line = filePDB.addEndLine(system.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom[v_num_atom], v_serial_number)
                v_file.write(line)
            v_file.close()
            print("file " + fileName + " created" + '\n')
######################################################################################################################
######################################################################################################################
    @staticmethod
    def addLine(atom, serialnumber,language):

        """Add line in PDB file with atom values.

        Input-----------------------
        atom (atom): atom object.
        serialnumber (int): value for atom's serial number.
        language (str): language to print.
        Output----------------------
        line (str): Line with atom info in PDB format.
        """
        num = serialnumber   
        if (abs(atom.x) > 999.99) or (abs(atom.x) > 999.99) or(abs(atom.x) > 999.99):
            print("Format error in real number " + atom.name + ' ' + '{:5d}'.format(serialnumber))
            error = messages.get_text_message('026', language)
            error = string.replace(error, '&', '{:5d}'.format(serialnumber), 1)
            raise IOError(error)
        if (serialnumber>99999):
            num = np.remainder(serialnumber,100000)

        line = '{:6}'.format(atom.tipus) + \
               '{:5d}'.format(num) + ' ' + \
               '{:^4}'.format(atom.name) + ' ' + \
               '{:3}'.format(atom.residue_name) + ' ' + \
               '{:1}'.format(atom.chain_id) + \
               '{:4}'.format(atom.residue_number) + '    ' + \
               '{: 8.3f}'.format(atom.x) + \
               '{: 8.3f}'.format(atom.y) + \
               '{: 8.3f}'.format(atom.z) + \
               '{: 6.2f}'.format(atom.occupancy) + \
               '{:6.2f}'.format(atom.temperature) + '          ' + \
               '{:2s}'.format(atom.element_symbol) + \
               '{:2s}'.format(atom.charge) + '\n'
        return line
######################################################################################################################
######################################################################################################################
    @staticmethod
    def addEndLine(atom, serialnumber):
        """Adds end line in PDB file with atom values.
        Input-----------------------
        atom (atom): atom object.
        serialnumber (int): value for atom's serial number.
        Output----------------------
        line (str): Line with TER in PDB format.
        """
        
        line = '{:6}'.format("TER") + \
               '{:5d}'.format(serialnumber) + '      ' + \
               '{:3}'.format(atom.residue_name) + ' ' + \
               '{:1}'.format(atom.chain_id) + \
               '{:4}'.format(atom.residue_number) + '\n'
        return line
######################################################################################################################
######################################################################################################################
    @staticmethod
    def addEndLine2():
        """Adds a plain line in PDB file with atom values.

        Output----------------------
        line (str): Line with TER in PDB format.
        """
        line = '{:6}'.format("TER") + '\n'
        return line
######################################################################################################################
