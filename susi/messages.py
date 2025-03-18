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
# without specific, prior, written permission. º
#########################################################################

class messages:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    @staticmethod
    def get_text_copyright():
        """Print the copyright information

        Output----------------------
        (str): copyright message.
        """

        msg = """ 
######################
# -----  SuSi   ----- 
# Structure Simulation
#  Copyright (C) 2025 
######################

"""
        return msg

######################################################################################################################
######################################################################################################################
    @staticmethod
    def get_text_license():
        """Print the license information

        Output----------------------
        (str): license message.
        """

        msg = """ 
######################################################################## 
#
#                 -----  SuSi   ----- 
#                 Structure Simulation
#                 Copyright (C) 2025 
#
#          David Naranjo        david.alejandro.naranjo@upc.edu 
#          Carlos Alemán        carlos.aleman@upc.edu
#          Juan Torras          joan.torras@upc.edu 
#
######################################################################## 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#        
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#          
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
########################################################################

"""
        return msg

######################################################################################################################
######################################################################################################################
    @staticmethod
    def get_text_message(numMessage, language):
        """Get the list of errors.

        Input-----------------------
           numMessage (str): Message key number.
           language (str): Language for printing.   
        Output-----------------------
           options[numMessage] (str): Error message.
        """
        if language == 'EN':
            options = {'001' : 'The input file can''t have less than 8 lines.',
                      '002' : 'Collisions could not be avoided for this molecule.',
                      '003' : 'system_name not found in the input file.',
                      '004' : 'box dimension not found in the input file.',
                      '005' : 'periodical_coordinates not found in the input file.',
                      '006' : 'minimum_distance not found in the input file.',
                      '007' : 'retries not found in the input file.',
                      '008' : 'maximum_deviation not found in the input file.',
                      '009' : 'molecules_number not found in the input file.',
                      '010' : 'PATH of prepi files not found in the input file.',
                      '011' : 'Error in X dimension units.',
                      '012' : 'Error in Y dimension units.',
                      '013' : '',
                      '014' : 'Error in Z dimension units.',
                      '015' : 'Error in minimum distance units.',
                      '016' : 'System name not found in file.',
                      '017' : 'Error in value of periodical coordinates. Possible values are Y or N.',
                      '018' : 'Error in minimum distance.',
                      '019' : 'Molecules number not found in the input file.',
                      '020' : 'PATH of prepi files not found in the input file.',
                      '021' : 'Molecules number in the input file is different from molecules number found.',
                      '022' : 'Error in box.',
                      '023' : 'Could not read file: &',
                      '024' : 'No errors found in the molecules building.',
                      '025' : 'Could not create file: &',
                      '026' : 'Large real number to fit a PDB format in atom &',
                      '027' : 'Wrong number of starting and closing atoms for crosslinking.',
                      '028' : 'Starting or closing atom does not exist in PREPI files. Wrong name for atom: &',
                      '029' : 'Error in distortion angle. Wrong or not specified value.', 
                      '030' : 'crosslinks_number not found in the input file.',
                      '031' : 'Crosslinks number in the input file is different from crosslinks number found.',
                      '032' : "Could not create cross links. The linear chains have been written in the restart file.",
                      '033' : 'System created successfully.',
                      '034' : 'Error in crosslink starting bond lenght. Wrong or not specified value.',
                      '035' : 'Error in the units. Only Angstrom (A) are supported.',
                      '036' : 'Minimum distance is greater than crosslink starting bond lenght.',
                      '037' : 'Warning: Minimum distance is greater than that of a typical C-C single bond.', 
                      '038' : 'Wrong number of removable atoms for starting and/or closing atoms for crosslinking. If no atoms need to be removed, specify as null.',
                      '039' : 'Removable atoms for the starting or closing atom does not exist in PREPI files. Wrong name for atom: &',
                      '040' : 'Removable atoms not bonded to the starting or closing atom. Wrong name for atom: &',
                      '041' : 'Error in crosslink closing bond lenght. Wrong or not specified value.', 
                      '042' : 'Minimum distance is greater than crosslink closing bond lenght.',
                      '043' : 'No closing atoms were found within the distance defined by the crosslink size.',
                      '044' : '"graft" can only be set as a closing atom.',
                      '045' : '"dendron" can only be set in a removable closing atom position.',
                      '046' : '"dendron" can only be used with "graft".'}
            return options[numMessage]
