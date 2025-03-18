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

import string 

class molecule:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self,name,ID):
        """Class representig molecules formed by residues objects.  
        
        Input-----------------------
            name (str): Molecule name.
            ID (str): Molecule ID.
        """
        self.t_residues = []
        self.name       = name
        self.ID         = ID
######################################################################################################################
######################################################################################################################
    def addResidue(self, residue):
        """Add a new residue objects in the residue list.  
        
        Input-----------------------
            residue (residue): Residue object.
        """        
        self.t_residues.append(residue)
######################################################################################################################
######################################################################################################################
    def printMolecule(self):
        """Prints all the residues of a molecule.
        """
        print( 'MOLECULE: Name:', self.name, ' ID:', self.ID, "Number of residues: ", str(len(self.t_residues)))
        for v_num_residue in range(0,len(self.t_residues)):
            self.t_residues[v_num_residue].printResidue()
##################################################################################################################
######################################################################################################################
    @staticmethod
    def getNextChainID(currentChainID):
        """Get next value in alphabet for a chain.
        
        Input-----------------------
        currentChainID (str): letter representing the name of the chain.
        Output----------------------
        nextChainID (str): letter representing the name of the chain.
        """

        if currentChainID == '' or currentChainID == 'Z': nextChainID = 'A'
        else:
            t_alphabet = list(string.ascii_uppercase)
            v_value = t_alphabet.index(currentChainID)
            nextChainID = t_alphabet[v_value+1]
        return nextChainID
######################################################################################################################
