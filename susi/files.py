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

import os.path   # To get file extensions with os.path.splitext
from .filePDB   import filePDB
from .filePREPI import filePREPI
######################################################################################################################
class files:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self,fileName):
        """File class constructor

        Input-----------------------
        fileName (str): Name of current file.
        """
        self.fileName = fileName
#Get file extension
        self.fileType = os.path.splitext(fileName)[1]
        if self.fileType in [".pdb", ".PDB"]:
            self.filePDB = filePDB()
        else:
            if self.fileType in [".prepi", ".PREPI"]:
                self.filePREPI = filePREPI()
            else:
                print("Error in format file")
######################################################################################################################
######################################################################################################################
    def readFile(self,systemName = "only for PDB"):
        """Read current file.

        Input-----------------------
        systemName (str): Name of current system.
        """
        if self.fileType in [".pdb", ".PDB"]:
            v_system = self.filePDB.readFile(self.fileName,systemName)
            return v_system
        else:
            if self.fileType in [".prepi", ".PREPI"]:
                self.filePREPI.readFile(self.fileName)
# Returns a list with PREPI lines
                v_PREPIlines = self.filePREPI.t_file_lines
                return v_PREPIlines
            else:
                print("Error in format file")
######################################################################################################################
######################################################################################################################
    @staticmethod
    def storeFile(fileName, system, language):
        """Stores current file.

        Input-----------------------
        fileName (str): Name of current file.
        system (system): System object.
        language (str): Language to print.
        """

        v_fileType = os.path.splitext(fileName)[1]
        if v_fileType in [".pdb", ".PDB"]:
            filePDB.storeFile(fileName,system,language)
        else:
            if v_fileType in [".prepi", ".PREPI"]:
                filePREPI.storeFile()
            else:
                print("Error in format file")
######################################################################################################################
######################################################################################################################
    def printFile(self):
        """Prints current file's info.
        """
        v_num_lines = 0
        if self.fileType in [".pdb", ".PDB"]:
            while v_num_lines < len(self.filePDB.t_file_lines ):
                print(self.filePDB.t_file_lines[v_num_lines])
                v_num_lines +=1
        else:
            if self.fileType in [".prepi", ".PREPI"]:
                while v_num_lines < len(self.filePREPI.t_file_lines ):
                    print(self.filePREPI.t_file_lines[v_num_lines])
                    v_num_lines +=1
            else:
                print("Error in format file")
######################################################################################################################
