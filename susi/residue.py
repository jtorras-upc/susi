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

class residue:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self, t_atom, name = "dummy", number = "dummy", moleculeID = "dummy"):
        """Class representig residue formed by atom objects.  
        
        Input-----------------------
            t_atom (list): List with atom objects.
            name (str): Residue name.
            number (str): Residue number.
            moleculeID (str): Molecule ID.
        """
        self.t_atom     = t_atom
        self.name       = name
        self.number     = number
        self.moleculeID = moleculeID
######################################################################################################################
######################################################################################################################
    def printResidue(self):
        """Prints all atoms of this residue.
        """

        print('RESIDUE NAME:', self.name, 'number:', self.number, ' molecule ID:', self.moleculeID)
        print('ATOMS: Number:' + str(len(self.t_atom)))
        for v_num_atom in range(0,len(self.t_atom)):
            self.t_atom[v_num_atom].printAtom()
            v_num_atom += 1
######################################################################################################################
######################################################################################################################
    def addAtom(self,atom):
        """Adds a new atom to current residue.

        Input-----------------------
            atom (atom): atom object.        
        """
        self.t_atom.append(atom)
######################################################################################################################
######################################################################################################################
    def getAtom(self,atom_name):
        """Gets atom from atribute t_atom of current residue

        Input-----------------------
            atom_name (str): Atom name. 
        """

        v_num_atom = 0
        while v_num_atom < len(self.t_atom):
            if self.t_atom[v_num_atom].name == atom_name:
                return self.t_atom[v_num_atom]
            v_num_atom += 1
        print("Atom " + atom_name + " not found")
######################################################################################################################
