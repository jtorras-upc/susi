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

class system:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################

    def __init__(self,name, minimum_distance, maximum_deviation, max_retries, box, language,state):
        """Class of objects that contains all the info of the system to build.

        Input-----------------------
            name (str): System name.
            minimum_distance (list): Size and unit of minimum distance.
            maximum_deviation (float): Value for maximum deviation.
            max_retries (int): Maximum number of retries to find a residue angle that don't collide with other atoms.
            box (list): List with 3 box dimensions (X,Y,Z) and unit[0] X dimension. [1] X unit [2] Y dimension. [3] Y unit [4] Z dimension. [5] Z unit.
            language (str): Language to print.
        """

        self.t_molecules       = []
        self.name              = name
        self.minimum_distance  = minimum_distance
        self.maximum_deviation = maximum_deviation
        self.max_retries       = max_retries
        self.box               = box
        self.language          = language
        self.state             = state
######################################################################################################################
######################################################################################################################
    def addMolecule(self,molecule):
        """Adds a new molecule to the system.

        Input-----------------------
            molecule (molecule): molecule object.        
        """
        self.t_molecules.append(molecule)
######################################################################################################################
######################################################################################################################
    def printSystemNumberElements(self):
        """Print the info of the current system.
        """             
        for v_num_record in range(0,len(self.t_molecules)):
            print('Molecule ID:', self.t_molecules[v_num_record].ID)
            if len(self.t_molecules[v_num_record].t_residues) > 0:
                print('Num Residues: %d' %len(self.t_molecules[v_num_record].t_residues))
            else:
                print('Num Residues: 0')
            print('Num Atoms: %d'    %self.t_molecules[v_num_record].getNumAtoms())
######################################################################################################################
######################################################################################################################
    def printSystemMolecules(self):
        """Print molecules of current system
        """
        print("SYSTEM NAME  :" + self.name + " Molecules number: " + str(len(self.t_molecules)))
        for v_num_molecule in range(0,len(self.t_molecules)):
            self.t_molecules[v_num_molecule].printMolecule()
######################################################################################################################
