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

class atom:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self,tipus,serial_number,name,chain_id, x,y,z,occupancy,temperature,element_symbol,charge, \
                 residue_name,residue_number,topology_type,internal_numbering,connectivity):
        """Class that represents the atoms objects for PDB and PREP formats.
        
        Input-----------------------
        tipus (str): Atom tipus.
        serial_number (int): Atom serial number.
        name (str): Atom name.
        chain_id (str): Value for chain id
        posX (float): X-axis coordinate.
        posY (float): Y-axis coordinate.
        posZ (float): Z-axis coordinate.
        occupancy (float): Value for occupancy.
        temperature (float): Value for temperature
        element_symbol (str): Value for element symbol.
        charge (float): Value for charge.
        residue_name (str): residue name.
        residue_number (int): value for residue number.
        """
        
        self.tipus                   = tipus
        self.serial_number           = serial_number
        self.name                    = name
        self.chain_id                = chain_id
        self.x                       = float(x)
        self.y                       = float(y)
        self.z                       = float(z)
        self.occupancy               = float(occupancy)
        self.temperature             = float(temperature)
        self.element_symbol          = element_symbol
        self.charge                  = charge
        self.residue_name            = residue_name
        self.residue_number          = residue_number
        self.topology_type           = topology_type
        self.internal_numbering      = internal_numbering
        self.connectivity            = connectivity
######################################################################################################################
######################################################################################################################
    def printAtom(self):
        """Print the features of an atom object.
        """
        print( str(self.serial_number), ' ', str(self.name), ' ', self.chain_id, ' ', (self.x), ' ', \
              str(self.y), ' ', str(self.z), ' ', str(self.occupancy), ' ', str(self.temperature), ' ', \
              self.element_symbol, ' ',  str(self.charge))
######################################################################################################################
######################################################################################################################
    @staticmethod
    def deleteDUMMlines(t_atom):
        """Delete all DUMM lines from molecule list.
        Input-----------------------
         t_atom (list):List with atom objects.
        Output----------------------
        t_atom_without_dumm_lines (list):List with atom objects with no dummy lines.
        """
        t_atom_without_dumm_lines = []
        for i in range(0,len(t_atom)):
            if t_atom[i].name != 'DUMM': t_atom_without_dumm_lines.append(t_atom[i])
        return t_atom_without_dumm_lines
######################################################################################################################
