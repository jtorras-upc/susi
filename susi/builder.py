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

from .system       import system
from .molecule     import molecule
from .residue      import residue
from .atom         import atom
from .files        import files
from .int2car_M    import int2car
from .messages     import messages
from .collision    import collision
from .progressBar  import progressBar
from .parserSusi   import parserSusi
import sys          # Allow use of function sys.exit()
import os           # Allow use of function path.exists
import random       # Allow functions to get random values
import numpy as np  # Allow functions to calculate distance between 2 points
import copy
import math
from datetime import datetime
######################################################################################################################
class builder():
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self, imFileName, verb, imTest = '', restartFile=''):
        """System builder.
            Read file information
            1-Checks if the file exists, otherwise it will show an error.
            2-File must have more than 7 lines, because 7th line and next lines are where molecules are stored.
            3-Gets values stored in file with function getFileValues.
            4-Gets all residues from molecules lines. Delete repeated value and store information in t_allResidues list.
            5-Store in the list t_pathFilePrepi the PATH prepi+namefile.
            6-Reads all PREPI files stored in t_pathFilePrepi and create list t_filePREPI with filePREPI values.
            7-Create system with system name value from file.
            8-Read all molecules from file and process them one by one with function processMolecule.
            
        Input-----------------------
            imFilename (str): Name of the input file.
            imTest (str): 'X' test mode activated, file.pdb will not be created.' ' test mode deactivated, file.pdb will be created.
            verb (int): Verbosity level. 0 non Verbosity. >=1  Verbosity.
        """
        t_lines           = []          # List with lines from file. Each line is stored in list with values separated by ','
                                        # [[system_name,1ALA], [molecules_number,1], [path_prepi,prepis], [mol1,ALA]]
        box               = []          # Array for store box dimensions and unit
                                        # [0] X dimension. [1] X unit
                                        # [2] Y dimension. [3] Y unit
                                        # [4] Z dimension. [5] Z unit
        periodical_coordinates = ''
        minimumDistance        = []     # Minimum distance lineage. Position 0: Number, Position 1: unit. Sample [10,A] = 10 A
        retries                = ""     # Number of retries to find a residue angle that don't collide with other atoms
        maximum_deviation      = ""     # % maximum deviation
        systemName             = ""     # System name, value obtained later from file
        pathPrepi              = ""     # Path of PREPI files, value obtained later from file
        t_molecules            = []     # Molecule list obtained from file: Molecule name, res1, res2...
        t_errors               = []     # To store errors found reading file
        t_allResidues          = []     # List of all residue names
        t_pathFilePrepi        = []     # List with PATH+FILENAME of PREPI files
        t_filePREPI            = []     # List of PREPI files
        v_moleculeID           = ''     # Value is a sequential value beginning with A: Possible values are A, B, C, D, E, F, G, H, I, J,...Z
        t_log                  = []
        language               = 'EN'   # Sample: 'EN' = English
        self.verbosity         = 0      # Level of verbosity in the program execution (0 verbosity by default)
        self.logFile           = None   # File object to write all program execution comments 
        t_moleculeRetries      = []     # Number of retries of molecule
        t_moleculeID           = []
        t_crosslinkRetries     = []
        t_crosslinkID          = []
        printRestart           = False  # Print the restart file at the end of the execution
        bond_length_CL_start   = []      #Bond distance between the first atom of the crosslink and the starting atom [value,unit]
        bond_length_CL_close   = []      #Bond distance between the last atom of the crosslink and the closing atom [value,unit]


# 0-Create a Log file where to print out all message depending on the verbosity level
        try:
            if verb >= 1:
                self.verbosity = verb

            self.logFile = open("susi.log","w")
            self.printLogs(messages.get_text_license(),0)
        except IOError:
            t_errors.append(('025',"susi.log"))          #Could not create file: &
            builder.printErrors(t_errors, language)
            sys.exit()
 
# 1-Checks if the file exists, otherwise it will show an error
        try:
            v_file = open(imFileName, "r")
        except IOError:
            t_errors.append(('023', imFileName))         #Could not read file: &
            builder.printErrors(t_errors, language)
            sys.exit()

# 2-Gets all lines from file and store each line in t_lines list
        for v_line in v_file:
            t_lines.append(v_line)
        v_file.close()

# 3-Get values stored in file with the class parserSusi and its function getAndCheckFileValues
        p = parserSusi()
        box, periodical_coordinates, minimumDistance, retries, maximum_deviation, systemName, pathPrepi, \
                                        t_molecules, t_crosslinks, t_starting_closing_atoms, t_sub_starting_closing_atoms , distortion_angle, bond_length_CL_start, bond_length_CL_close, t_cycles, t_errors = p.getAndCheckFileValues(t_lines)
        #initialization of collision
        n_part = int(min(box[0],box[2],box[4])/(10*minimumDistance[0])) #Number of partitions of every axis for collision analysis 
        
        self.coll=collision(box,minimumDistance[0],n_part)    
        self.coll.createDic()
    
        if len(t_errors) > 0: # If errors found exit function
            builder.printErrors(t_errors, language)
            sys.exit()
        else:

# 4-Gets all residues from molecules lines and crosslinks lines. Delete repeated value and store information in t_allResidues list
            for v_number in range(0,len(t_molecules)):
                for v_number2 in range(0,len(t_molecules[v_number])):
                    t_allResidues.append(t_molecules[v_number][v_number2])
            for v_number in range(0,len(t_crosslinks)):
                for v_number2 in range(0,len(t_crosslinks[v_number])):
                    t_allResidues.append(t_crosslinks[v_number][v_number2])
            t_allResidues = sorted(list(dict.fromkeys(t_allResidues)))            # Order list and delete duplicates

# 5-Store in the list t_pathFilePrepi the PATH prepi+namefile: Example [prepis/ALA.prepi],[prepis/GLY.prepi]
            for v_number in range(0,len(t_allResidues)):
                v_pathFile = pathPrepi + '/' + t_allResidues[v_number] + '.prepi'
# Checks if prepi filename exists, otherwise store a error message in list t_errors
                try:
                    v_file = open(v_pathFile, "r")
                    t_pathFilePrepi.append(v_pathFile)
                    v_file.close()
                except IOError:
                    t_errors.append(('023', v_pathFile))         #Could not read file: &
# If errors are found, display in screen and leave function
            if len(t_errors) > 0:
                builder.printErrors(t_errors, language)
                sys.exit()
            else:
                
                random.seed(int(datetime.now().strftime("%Y%m%d%H%M%S")))
                #RANDOMSEED
                #random.seed(546)


# 6-Reads all PREPI files stored in t_pathFilePrepi and create list t_filePREPI with filePREPI values
                for v_number in range(0,len(t_pathFilePrepi)):
                    v_filePREPI = files(t_pathFilePrepi[v_number])
                    v_filePREPI.readFile()
                    t_filePREPI.append(v_filePREPI.filePREPI)
                #In case of crosslinks, check whether the names of the starting and closing atoms exist in the Prepi files
                dic_starting_closing = {}
                if len(t_starting_closing_atoms) != 0:
                    for pair_atoms in t_starting_closing_atoms:
                        control_name = 0
                        for name_atom in pair_atoms:                            
                            #control_name = 0
                            if name_atom == "graft":
                                control_name += 1 
                            else:     
                                dic_starting_closing[name_atom] = ""
                                t_list_name_atom = []
                                #control_name = 0
                                for file_prepi in t_filePREPI:
                                    if name_atom in file_prepi.t_name: 
                                        t_list_name_atom.append(file_prepi.name)        
                                        control_name += 1
                                t_list_name_atom = list(set(t_list_name_atom))
                                dic_starting_closing[name_atom] = t_list_name_atom
                            if control_name == 0: 
                                t_errors.append(('028', name_atom))         #Closing or starting atom does not exist in the prepi files
                                builder.printErrors(t_errors, language)
                                sys.exit()
                #In case of crosslinks, check whether the names of the removable starting and closing atoms exist in the Prepi files

                if len(t_sub_starting_closing_atoms) != 0:
                    for pair_atoms in t_sub_starting_closing_atoms:
                        for name_atom in pair_atoms:
                            #dic_starting_closing[name_atom] = ""
                            control_name = 0       
                            for file_prepi in t_filePREPI:
                                for name_prepi in file_prepi.t_name:
                                    if name_atom == "null" or name_atom == "branch" or name_atom == "dendron":
                                        control_name += 1
                                    elif name_atom == name_prepi: # name_atom in name_prepi: #if we need atom types
                                        control_name += 1   
                            if control_name == 0: 
                                t_errors.append(('039', name_atom))         #Closing or starting atom does not exist in the prepi files
                                builder.printErrors(t_errors, language)
                                sys.exit()

                #Check if the subtittution atoms are bonded to the starting or closing atom
                t_available_sub_starting_atom = []
                if len(t_sub_starting_closing_atoms) != 0:
                    for idx_pair_atom in range(0,len(t_sub_starting_closing_atoms)):
                        starting_atom = t_starting_closing_atoms[idx_pair_atom][0]
                        sub_starting_atom = t_sub_starting_closing_atoms[idx_pair_atom][0]
                        for file_prepi in t_filePREPI:                           
                            if starting_atom in file_prepi.t_name:
                                numbering_starting = file_prepi.t_numbering[file_prepi.t_name.index(starting_atom)]
                                t_sub_atom = []
                                for atom_sub in file_prepi.t_name:
                                    if sub_starting_atom != "null" or sub_starting_atom != "branch" :
                                        if sub_starting_atom == atom_sub:
                                            index_sub = file_prepi.t_name.index(atom_sub)
                                            connect_bond = file_prepi.t_con[index_sub][0] #number representing the atom to which the removable atom is bonded in the PREPI file
                                            if int(connect_bond) == int(numbering_starting):
                                                t_sub_atom.append(atom_sub)
                                if sub_starting_atom == "null":
                                    t_available_sub_starting_atom.append(["null"])    
                                elif sub_starting_atom == "branch":
                                    t_available_sub_starting_atom.append(["branch"])
                                else:
                                    t_available_sub_starting_atom.append(t_sub_atom)

                        for t_at in t_available_sub_starting_atom:
                            if len(t_at) == 0:
                                t_errors.append(('040',sub_starting_atom))         #removable atoms are not bonded to Closing or starting atom
                                builder.printErrors(t_errors, language)
                                sys.exit()

                t_available_sub_closing_atom = []
                if len(t_sub_starting_closing_atoms) != 0:
                    for idx_pair_atom in range(0,len(t_sub_starting_closing_atoms)):
                        closing_atom = t_starting_closing_atoms[idx_pair_atom][1]
                        sub_closing_atom = t_sub_starting_closing_atoms[idx_pair_atom][1]
                        for file_prepi in t_filePREPI:
                            if closing_atom in file_prepi.t_name:
                                numbering_starting = file_prepi.t_numbering[file_prepi.t_name.index(closing_atom)]
                                t_sub_atom = []
                                for atom_sub in file_prepi.t_name:
                                    if sub_closing_atom != "null" or sub_closing_atom != "branch":
                                        if sub_closing_atom == atom_sub:
                                            index_sub = file_prepi.t_name.index(atom_sub)
                                            connect_bond = file_prepi.t_con[index_sub][0] #number representing the atom to which the removable atom is bonded in the PREPI file
                                            if int(connect_bond) == int(numbering_starting):
                                                t_sub_atom.append(atom_sub)
                                if sub_closing_atom == "null":
                                    t_available_sub_closing_atom.append(["null"])    
                                elif sub_closing_atom == "branch":
                                    t_available_sub_closing_atom.append(["branch"])    
                                else:
                                    t_available_sub_closing_atom.append(t_sub_atom)

                        for t_at in t_available_sub_closing_atom:
                            if len(t_at) == 0: 
                                t_errors.append(('040',sub_closing_atom))         #removable atoms are not bonded to Closing or starting atom
                                builder.printErrors(t_errors, language)
                                sys.exit()

                #Create a list of [molecule, residue] indexes containg the starting or closing atoms
                dic_mol_res = {}
                for key_atom in dic_starting_closing:
                    t_mol_res_cl = []
                    for res_name in dic_starting_closing[key_atom]:
                        for idx_molecule in range (0, len(t_molecules)):
                            for idx_residue in range (0, len(t_molecules[idx_molecule])):
                                if t_molecules[idx_molecule][idx_residue] == res_name:
                                    t_mol_res_cl.append([idx_molecule,idx_residue])
                    dic_mol_res[key_atom] = t_mol_res_cl
                
                #Check if the graft keyword is only in a closing atom position

                for pair_atoms in t_starting_closing_atoms:
                    if pair_atoms[0] == "graft":
                        t_errors.append(('044'))         #graft keyword in starting position
                        builder.printErrors(t_errors, language)
                        sys.exit()
                               #Check if the graft keyword is only in a closing atom position
                
                #Check if the dendron keyword is only in a removable closing atom position
                for pair_atoms in t_sub_starting_closing_atoms:
                    if pair_atoms[0] == "dendron":
                        t_errors.append(('045'))         #dendron keyword in starting position
                        builder.printErrors(t_errors, language)
                        sys.exit()
                #Check if the dendron keyword is used with the graft keyword
                for ii in range (len(t_sub_starting_closing_atoms)):
                    if t_sub_starting_closing_atoms[ii][1] == "dendron" and  t_starting_closing_atoms[ii][1] != "graft":
                        t_errors.append(('046'))         #dendron keyword not used with graft
                        builder.printErrors(t_errors, language)
                        sys.exit()



# 7-Create system with system name value from file
                initial_state = [0,0,0]  #Number of molecules, number of resiudes, 1 if its a restart file 0 otherwise
                self.system_pbc  = system(systemName, minimumDistance, maximum_deviation, retries, box, language, initial_state)
                self.system_nopbc = system(systemName, minimumDistance, maximum_deviation, retries, box, language,initial_state)
                if restartFile:   
                    try:
                        self.system_nopbc=files(restartFile).readFile()
                        self.system_pbc=copy.deepcopy(self.system_nopbc)
                        for v_num_molecule in range (0,len(self.system_pbc.t_molecules)):   
                            for v_num_residue in range (0,len(self.system_pbc.t_molecules[v_num_molecule].t_residues)):   
                                for v_num_atom in range (0,len(self.system_pbc.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom)):
                                    a=self.system_pbc.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom[v_num_atom]  
                                    a.x, a.y, a.z = self.periodicBoundaryCondition(self.system_pbc.box,a.x,a.y,a.z) 
                                    self.system_pbc.t_molecules[v_num_molecule].t_residues[v_num_residue].t_atom[v_num_atom]=a  
                    except IOError:
                        t_errors.append(('023', restartFile))
                        builder.printErrors(t_errors, language)
                        sys.exit()
                
# 8-Read all molecules from file and process them one by one with function processMolecule
                for fillMoleculeRetries in range(0,len(t_molecules)): t_moleculeRetries.append(0)

                for moleculeNumber in range(0,len(t_molecules)):         #Fill t_moleculeID with next sequential value: A,B,C,D...)
                    v_moleculeID = molecule.getNextChainID(v_moleculeID)
                    t_moleculeID.append(v_moleculeID)

                # Create a progress bar of execution
                numResidues = 0
                for numMol in range(0,len(t_molecules)) : numResidues += len(t_molecules [numMol])
                for numMol in range(0,len(t_crosslinks)): numResidues += len(t_crosslinks[numMol])
                #width of the terminal screen
                width = int(os.get_terminal_size()[0])
                self.bar = progressBar(numResidues, "Computing","Residues", width-60) #60 is an arbitrary value for the size of the progress bar

                # Starting main loop through all system molecules
                moleculeNumber = self.system_nopbc.state[0] #Starts from 0
                self.system_pbc.state[0]=self.system_nopbc.state[0]
                v_exit = False
                total_molecules = len(t_molecules) #number of molecules not considering the crosslinks
                t_molecules_idx = [] # List with the numerical indexes of the molecules placed in the system
                              
                while moleculeNumber < total_molecules  and v_exit == False:
                    t_moleculeRetries[moleculeNumber] += 1
                    #Check if there are still retries for this molecule
                    if t_moleculeRetries[moleculeNumber] <= self.system_pbc.max_retries: 
                        
                        self.printLogs('Begin process for molecule '+ t_moleculeID[moleculeNumber],0)#0 is the verbosity level for printLogs
                        v_error = self.__addNewMolecule(t_molecules[moleculeNumber], moleculeNumber, \
                                                       t_moleculeID[moleculeNumber], t_filePREPI,t_molecules[moleculeNumber],t_cycles[moleculeNumber])
                        
                        #After __addNewMolecule. Molecule object is created. 
                        # If no collisions are detected resiude and atoms of that molecule are placed and created,and the output of the function is 0.
                        # If no residues are placed v_error is assigned with a tuple 
                        if v_error == 0:
                           t_molecules_idx.append(moleculeNumber)
                           moleculeNumber +=1               #Process the next molecule
                           self.system_nopbc.state[0] = moleculeNumber                    
                           self.system_pbc.state[0] = self.system_nopbc.state[0]

                           
                        else: #If resiudes could not be placed, we recover systems (pbc and no pbc) and 1
                           self.bar.substract(len(self.system_pbc.t_molecules[moleculeNumber].t_residues)) #Update Progress bar counter. #Delete as many elements as the number of residues of the current molecule
                           #If collisions are detected, no residues are placed, but the molecule object will remain.
                           # The molecule object with no residues is deleted and the retries increase. 
                           del self.system_pbc.t_molecules[-1]              #Delete molecule from system
                           del self.system_nopbc.t_molecules[-1]             #Delete molecule from system


                    else: #No retries reminding
                        #The flag v_exit will be activated. The molecules that were correctly placed in the system will never be deleted. 
                        t_errors.append('002')    #Collisions could not be avoided for this system
                        if moleculeNumber == 0: 
                            if imTest == '':
                                numTotalResidues=0
                                for v_num_molecule in range (0,len(self.system_nopbc.t_molecules)): #
                                    numTotalResidues +=len(self.system_nopbc.t_molecules[v_num_molecule].t_residues)
                                if numTotalResidues == 0:
                                    printRestart = False

                        else:
                            t_moleculeRetries[moleculeNumber] = 0    #Reset retries for this molecule 
                            moleculeNumber -= 1
                            self.system_nopbc.state[0] = moleculeNumber                    
                            self.system_pbc.state[0] = self.system_nopbc.state[0]
                            self.system_nopbc.state[2] = 1  #The reset flag will be placed in the pdb file
                            self.system_pbc.state[2] = 1  #The reset flag will be placed in the pdb file
                            printRestart = True

                        v_exit = True
# 9-Process and creates the crosslinks if exist
            if len(t_crosslinks) and not v_exit > 0:
                #Create a file with the information of the new bonds between the molecules and the crossllinks
                crosslink_info = open(self.system_pbc.name + '_CLbond.dat','w')
#####################################################################################################
                if len(t_errors) or v_exit > 0: 
                    self.endSusi(t_errors,imTest,periodical_coordinates,language,box,minimumDistance, \
                         retries,maximum_deviation,pathPrepi,t_molecules,t_allResidues,t_pathFilePrepi,printRestart)

                for fillCrosslinksRetries in range(0,len(t_crosslinks)): t_crosslinkRetries.append(0)  #Start loop from crosslinking
                for crosslinksNumber in range(0,len(t_crosslinks)):         #Fill t_crosslinkID with next sequential value: A,B,C..
                    v_moleculeID = molecule.getNextChainID(v_moleculeID)
                    t_crosslinkID.append(v_moleculeID)
                #Starting the loop for the crosslink
                v_exit = False 

                t_dis_cl = []
                t_v_growth_cl = []
                t_t_residues_cl = []
                
                #Measure the size of all the crosslink molecules
                for cl_idx in range (len(t_crosslinks)):
                    dis , v_growth , t_residues = self.__MeasureMolecule(t_crosslinks[cl_idx],len(t_molecules)+cl_idx,t_crosslinkID[cl_idx], t_filePREPI)
                    t_dis_cl.append(dis)
                    t_v_growth_cl.append(v_growth)
                    t_t_residues_cl.append(t_residues)
            
                #Fill retries for pair [molecule, residue]
                dic_mol_res_retries = {}
                for key_mol_res in dic_mol_res:
                    t_mol_res_retries = []
                    for _ in range(0,len(dic_mol_res[key_mol_res])): 
                        t_mol_res_retries.append(0)
                    dic_mol_res_retries[key_mol_res] = t_mol_res_retries   
                
                crosslinkNumber = 0
                #Building the crosslinks within the system                          
                while  crosslinkNumber < len(t_crosslinks) and v_exit == False :
                    t_crosslinkRetries[crosslinkNumber] += 1
                    starting_atom_name = t_starting_closing_atoms[crosslinkNumber][0]
                    closing_atom_name  = t_starting_closing_atoms[crosslinkNumber][1]
                    dis_cl = t_dis_cl[crosslinkNumber]
                    v_growth_cl = t_v_growth_cl[crosslinkNumber]
                    t_residues_cl = t_t_residues_cl[crosslinkNumber]
                    
                    if closing_atom_name != "graft": 
                        tolerance = 1000 
                    else: 
                        tolerance =1

                    if t_crosslinkRetries[crosslinkNumber] <= self.system_pbc.max_retries*tolerance:
                        totalResidues_cl = len (t_residues_cl)
                        
                        if t_crosslinkRetries[crosslinkNumber] == 1 and closing_atom_name != "graft":
                            t_possibilities = []
                            for element in dic_mol_res[starting_atom_name]:
                                rnd_mol_ini = element[0] 
                                rnd_res_ini = element[1]                            
                                HM_atom_ini = self.__get_atom_ini_idx(rnd_mol_ini,rnd_res_ini,starting_atom_name)#Number representing the index of an starting atom
                                if HM_atom_ini != None:
                                    #atom_ini_sub = self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[HM_atom_ini] #atom object (Starting atom)
                                    atom_ini_sub = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[HM_atom_ini] #atom object (Starting atom)
                                    
                                    #Seleccting the closing atom located at a distance determined by the length of crosslinker
                                    mol_end, res_end, HM_atom_end =  self.__get_atom_end_idx( dic_mol_res[closing_atom_name], dis_cl, bond_length_CL_start[0],bond_length_CL_close[0], atom_ini_sub,rnd_mol_ini, closing_atom_name, False)
                                    if HM_atom_end != None:
                                        t_possibilities.append([rnd_mol_ini,rnd_res_ini,HM_atom_ini,mol_end, res_end, HM_atom_end])
                        
                        #########################################START BUILDING GRAFTS/DENDRIMERS################################################        
                        if closing_atom_name == "graft":
                            if  t_sub_starting_closing_atoms[crosslinkNumber][1] != "dendron" or (t_sub_starting_closing_atoms[crosslinkNumber][1] == "dendron" and t_sub_starting_closing_atoms[crosslinkNumber][0] == "null"):
                                rnd_mol_ini,rnd_res_ini = random.choice(dic_mol_res[starting_atom_name])
                                HM_atom_ini = self.__get_atom_ini_idx(rnd_mol_ini,rnd_res_ini,starting_atom_name)#Number representing the index of a starting atom
                            else:
                                it_dendron_res = []
                                for mol_res in dic_mol_res[starting_atom_name]:
                                    for atm_rmv in self.system_pbc.t_molecules[mol_res[0]].t_residues[mol_res[1]].t_atom:
                                        if atm_rmv.name == t_sub_starting_closing_atoms[crosslinkNumber][0]:
                                            it_dendron_res.append(mol_res)
                                if len (it_dendron_res) != 0:
                                    rnd_mol_ini,rnd_res_ini = random.choice(it_dendron_res) 
                                    HM_atom_ini = self.__get_atom_ini_idx(rnd_mol_ini,rnd_res_ini,starting_atom_name)#Number representing the index of a starting atom
                            
                            #extract the growth vector where the starting atom is located
                            atm_idx = []
                            for atom_elm in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                atm_idx.append(atom_elm.topology_type)
                            atom_M0_idx_start = atm_idx.index('M')
                            atom_Mlast_idx_start = len(atm_idx) - 1 - atm_idx[::-1].index('M')   

                            atom_M0_start = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[atom_M0_idx_start]
                            atom_Mlast_start = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[atom_Mlast_idx_start]

                            v_growth_res_start = np.array([atom_Mlast_start.x,atom_Mlast_start.y,atom_Mlast_start.z]) - np.array([atom_M0_start.x,atom_M0_start.y,atom_M0_start.z])
                            v_growth_res_start /= np.linalg.norm(v_growth_res_start)


                            ####Removing from collision dictionary the heavy atoms when its specified
                            if t_sub_starting_closing_atoms[crosslinkNumber][0] != "null" and t_sub_starting_closing_atoms[crosslinkNumber][0] != "branch":    
                                starting_idx = 0
                                for atom_rmv in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                    if atom_rmv.name == t_sub_starting_closing_atoms[crosslinkNumber][0]:
                                        atom_to_remove_start = atom_rmv #atom obect to be subsituted in the starti position (sub start atom)
                                        removable_starting_idx = starting_idx
                                        self.coll.removeFromDictionary([atom_to_remove_start])
                                    starting_idx += 1
                            elif t_sub_starting_closing_atoms[crosslinkNumber][0] == "branch":
                                #Atoms list of the starting atom residue. Does not contain DUMMY atoms
                                t_atoms_starting = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom
                                #Get indexes of atoms connected to the starting atom in the residue
                                t_connected_starting, in_starting_atom = self.__find_connected_atoms(t_atoms_starting,starting_atom_name)
                                #Get a list of non-M heavy atoms connected to the starting atom
                                t_removable_starting = []
                                t_removable_starting_idx = []
                                for i in t_connected_starting:
                                    if t_atoms_starting[i].topology_type != 'M':
                                        if 'H' not in t_atoms_starting[i].name:
                                            t_removable_starting.append(t_atoms_starting[i])
                                            t_removable_starting_idx.append(i)
                                        elif t_atoms_starting[i].connectivity[0] != in_starting_atom:
                                            t_removable_starting.append(t_atoms_starting[i])
                                            t_removable_starting_idx.append(i)                    
                                #Removing from the collision dictionary the heavy non-M atoms connected to the starting atom
                                self.coll.removeFromDictionary(t_removable_starting)
                                t_removable_starting_idx.sort(reverse=True)

                            atom_ini_sub = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[HM_atom_ini] #atom object (starting atom)
                            
                            v_error = self.__addNewGraft(t_crosslinks[crosslinkNumber], len(t_molecules)+crosslinkNumber,t_crosslinkID[crosslinkNumber], t_filePREPI,rnd_mol_ini, rnd_res_ini,HM_atom_ini,v_growth_res_start,v_growth_cl,atom_ini_sub,t_residues_cl,bond_length_CL_start[0])

                            if v_error == 0:
                                
                                #Get the indexes of the residues in LEaP format
                                if  t_sub_starting_closing_atoms[crosslinkNumber][1] == "dendron": 
                                    dendron_flag = True
                                else:
                                    dendron_flag = False

                                res_ini_idx_tleap, res_end_idx_tleap, CL_ini_idx_tleap, CL_end_idx_tleap = self.__getLEaPindex(rnd_mol_ini, rnd_res_ini,t_molecules,None,None,crosslinkNumber,t_crosslinks, dendron_flag)

                                atom_ini = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[HM_atom_ini]
                                
                                CL_atom_ini = self.system_pbc.t_molecules[len(t_molecules) + crosslinkNumber].t_residues[0].t_atom[0]
                            
                                crosslink_info.write("#####  Bond info GRAFT: " + str(crosslinkNumber + 1) + "#### \n")    
                                
                                #Delete the removable atoms if specified for starting atom
                                if t_sub_starting_closing_atoms[crosslinkNumber][0] == "branch":
                                    for idx_rm_start in t_removable_starting_idx:
                                        del self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[idx_rm_start]
                                        del self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[idx_rm_start]
                                    
                                    #Rename the residue when atoms have been removed
                                    old_res_name = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].name
                                    new_res_name = old_res_name[0:2] + "b" #Format 2 letters + b 

                                    for atom_new in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name             
                                    for atom_new in self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name                 

                                    #######################################################################################################################
                                    crosslink_info.write("#Removed Branch Atoms. Renamed residue: " +  old_res_name + ".    New residue name: " + new_res_name + "\n")
                                    #######################################################################################################################                                                                                                                   

                                elif t_sub_starting_closing_atoms[crosslinkNumber][0] != "null":
                                    #Remove the removable atom from the starting residue

                                    del self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[removable_starting_idx]
                                    del self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[removable_starting_idx]

                                    #Rename the residue when atoms have been removed
                                    old_res_name = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].name
                                    new_res_name = old_res_name[0:2] + atom_to_remove_start.name[-1]

                                    for atom_new in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name                                    
                                    for atom_new in self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name                                    
                                    
                                    #######################################################################################################################
                                    crosslink_info.write("#Removed Atom: " + atom_to_remove_start.name + ". Renamed residue: " +  old_res_name +  ".    New residue name: " + new_res_name + "\n")
                                    #######################################################################################################################                                                                                            
                                
                                crosslink_info.write("bond    mol." + str(res_ini_idx_tleap) + "." + atom_ini.name + "    mol." + str(CL_ini_idx_tleap) + "." + CL_atom_ini.name + "\n")                                
                                crosslink_info.write("\n")

                                #Process the next crosslinker molecule
                                if  t_sub_starting_closing_atoms[crosslinkNumber][1] == "dendron":
                                    
                                    new_chain_idx = len(t_molecules)+crosslinkNumber
                                    idx_residue_CL = 0
                                    for new_residue in self.system_pbc.t_molecules[new_chain_idx].t_residues:
                                        for new_atom in new_residue.t_atom: 
                                            dic_mol_res.setdefault(new_atom.name,[]).append([new_chain_idx,idx_residue_CL])
                                        idx_residue_CL += 1


                                if  t_sub_starting_closing_atoms[crosslinkNumber][1] != "dendron":
                                    #Remove the molecule and the residue of starting and closing atoms of the [molecule,residue] so that they are not consider again for new crosslinks
                                    dic_mol_res[starting_atom_name].remove([rnd_mol_ini, rnd_res_ini])
                
                                crosslinkNumber += 1        
                                self.system_nopbc.state[0] = len(t_molecules) + crosslinkNumber                    
                                self.system_pbc.state[0] = self.system_nopbc.state[0]
                                ##################################################
                                    
                            else:
                                self.bar.substract(len(self.system_pbc.t_molecules[len(t_molecules)+crosslinkNumber].t_residues))
                                del self.system_pbc.t_molecules[-1]
                                del self.system_nopbc.t_molecules[-1]
                                
                                if t_sub_starting_closing_atoms[crosslinkNumber][0] == "branch":
                                    for elem in t_removable_starting:
                                        self.coll.addAtomToDictionary(elem)
                                
                                elif t_sub_starting_closing_atoms[crosslinkNumber][0] != "null":
                                        self.coll.addAtomToDictionary(atom_to_remove_start)
                        #########################################END BUILDING GRAFTS/DENDRIMERS##################################################        
                        
                        #########################################START BUILDING CROSSLINKS#######################################################        

                        elif closing_atom_name != "graft" and len(t_possibilities) != 0:
                            ##################################
                            pair_possible =random.choice(t_possibilities)
                            rnd_mol_ini = pair_possible[0]
                            rnd_res_ini = pair_possible[1]
                            HM_atom_ini = pair_possible[2]
                            mol_end = pair_possible[3]
                            res_end = pair_possible[4]
                            HM_atom_end = pair_possible[5]

                            index_mol_res = dic_mol_res[starting_atom_name].index([rnd_mol_ini, rnd_res_ini ])
                            dic_mol_res_retries[starting_atom_name][index_mol_res] += 1 

                            atom_ini_sub = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[HM_atom_ini] #atom object (starting atom)

                            atom_end_sub = self.system_pbc.t_molecules[mol_end].t_residues[res_end].t_atom[HM_atom_end] #atom object (Closing atom)

                            atom_end_idx = [mol_end, res_end, HM_atom_end]
                            #Compute the vector direction between the initial and final atoms
                            v_end_ini = np.array((atom_end_sub.x, atom_end_sub.y, atom_end_sub.z))-np.array((atom_ini_sub.x, atom_ini_sub.y, atom_ini_sub.z))
                            #Transform in unitary vector
                            v_end_ini_sub = v_end_ini/np.linalg.norm(v_end_ini)

                            ####Removing from collision dictionary the heavy atoms when its specified
                            if t_sub_starting_closing_atoms[crosslinkNumber][0] != "null" and t_sub_starting_closing_atoms[crosslinkNumber][0] != "branch":    
                                starting_idx = 0
                                for atom_rmv in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                    if atom_rmv.name == t_sub_starting_closing_atoms[crosslinkNumber][0]:
                                        atom_to_remove_start = atom_rmv #atom obect to be subsituted in the starti position (sub start atom)
                                        removable_starting_idx = starting_idx
                                        self.coll.removeFromDictionary([atom_to_remove_start])
                                    starting_idx += 1
                            elif t_sub_starting_closing_atoms[crosslinkNumber][0] == "branch":
                                #Atoms list of the starting atom residue. Does not contain DUMMY atoms
                                t_atoms_starting = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom
                                #Get indexes of atoms connected to the starting atom in the residue
                                t_connected_starting, in_starting_atom = self.__find_connected_atoms(t_atoms_starting,starting_atom_name)
                                #Get a list of non-M heavy atoms connected to the starting atom
                                t_removable_starting = []
                                t_removable_starting_idx = []
                                for i in t_connected_starting:
                                    if t_atoms_starting[i].topology_type != 'M':
                                        if 'H' not in t_atoms_starting[i].name:
                                            t_removable_starting.append(t_atoms_starting[i])
                                            t_removable_starting_idx.append(i)
                                        elif t_atoms_starting[i].connectivity[0] != in_starting_atom:
                                            t_removable_starting.append(t_atoms_starting[i])
                                            t_removable_starting_idx.append(i)                    
                                #Removing from the collision dictionary the heavy non-M atoms connected to the starting atom
                                self.coll.removeFromDictionary(t_removable_starting)
                                t_removable_starting_idx.sort(reverse=True)
                                                       
                            v_error, new_mol_end, new_res_end, new_HM_atom_end = self.__addNewCrossLink(t_crosslinks[crosslinkNumber], len(t_molecules)+crosslinkNumber, \
                                                                                                    t_crosslinkID[crosslinkNumber], t_filePREPI, v_growth_cl, v_end_ini_sub, \
                                                                                                        atom_ini_sub,atom_end_idx, closing_atom_name,t_residues_cl, \
                                                                                                            totalResidues_cl,dic_mol_res[closing_atom_name],rnd_mol_ini, \
                                                                                                                distortion_angle, bond_length_CL_start[0], bond_length_CL_close[0],t_sub_starting_closing_atoms[crosslinkNumber][1])
                            
                            if v_error == 0:
                                
                                #Get the indexes of the residues in LEaP format

                                res_ini_idx_tleap, res_end_idx_tleap, CL_ini_idx_tleap, CL_end_idx_tleap = self.__getLEaPindex(rnd_mol_ini, rnd_res_ini,t_molecules,new_mol_end,new_res_end,crosslinkNumber,t_crosslinks,False)

                                atom_ini = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[HM_atom_ini]
                                
                                CL_atom_ini = self.system_pbc.t_molecules[len(t_molecules) + crosslinkNumber].t_residues[0].t_atom[0]

                                atom_end = self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[new_HM_atom_end]
                                
                                residue_end_CL = self.system_pbc.t_molecules[len(t_molecules) + crosslinkNumber].t_residues[-1]
                                
                                CL_atom_end  = residue_end_CL.t_atom[self.__getm1atomResidue(residue_end_CL,t_filePREPI)-3]
                                
                                crosslink_info.write("#####  Bond info CL: " + str(crosslinkNumber + 1) + "#### \n")    
                                
                                #Delete the removable atoms if specified for starting atom
                                if t_sub_starting_closing_atoms[crosslinkNumber][0] == "branch":
                                    for idx_rm_start in t_removable_starting_idx:
                                        del self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[idx_rm_start]
                                        del self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[idx_rm_start]
                                    
                                    #Rename the residue when atoms have been removed
                                    old_res_name = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].name
                                    new_res_name = old_res_name[0:2] + "b" #Format 2 letters + b 

                                    for atom_new in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name             
                                    for atom_new in self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name                 

                                    #######################################################################################################################
                                    crosslink_info.write("#Removed Branch Atoms. Renamed residue: " +  old_res_name + ".    New residue name: " + new_res_name + "\n")
                                    #######################################################################################################################                                                                                                                   

                                elif t_sub_starting_closing_atoms[crosslinkNumber][0] != "null":
                                    #Remove the removable atom from the starting residue
                                    del self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[removable_starting_idx]
                                    del self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom[removable_starting_idx]

                                    #Rename the residue when atoms have been removed
                                    old_res_name = self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].name
                                    new_res_name = old_res_name[0:2] + atom_to_remove_start.name[-1]

                                    for atom_new in self.system_pbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name                                    
                                    for atom_new in self.system_nopbc.t_molecules[rnd_mol_ini].t_residues[rnd_res_ini].t_atom:
                                        atom_new.residue_name = new_res_name                                    
                                    
                                    #######################################################################################################################
                                    crosslink_info.write("#Removed Atom: " + atom_to_remove_start.name + ". Renamed residue: " +  old_res_name +  ".    New residue name: " + new_res_name + "\n")
                                    #######################################################################################################################                                                                                            
                                
                                crosslink_info.write("bond    mol." + str(res_ini_idx_tleap) + "." + atom_ini.name + "    mol." + str(CL_ini_idx_tleap) + "." + CL_atom_ini.name + "\n")
                                
                                #Delete the removable atoms if specified for closing atom
                                if t_sub_starting_closing_atoms[crosslinkNumber][1] == "branch":
                                    t_atoms_closing = self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom
                                    #Get indexes of atoms connected to the closing atom in the residue
                                    t_connected_closing, in_closing_atom = self.__find_connected_atoms(t_atoms_closing,closing_atom_name)
                                    #Get a list of non-M heavy atoms connected to the closing atom
                                    t_removable_closing = []
                                    t_removable_closing_idx = []
                                    for i in t_connected_closing:
                                        if t_atoms_closing[i].topology_type != 'M':
                                            if 'H' not in t_atoms_closing[i].name:
                                                t_removable_closing.append(t_atoms_closing[i])
                                                t_removable_closing_idx.append(i)
                                            elif t_atoms_closing[i].connectivity[0] != in_closing_atom:
                                                t_removable_closing.append(t_atoms_closing[i])
                                                t_removable_closing_idx.append(i)                                    
                                    t_removable_closing_idx.sort(reverse=True)
                                    for idx_rm_closing in t_removable_closing_idx:
                                        del self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[idx_rm_closing]                                    
                                        del self.system_nopbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[idx_rm_closing] 

                                    #Rename the residue when atoms have been removed
                                    old_res_name = self.system_nopbc.t_molecules[new_mol_end].t_residues[new_res_end].name
                                    new_res_name = old_res_name[0:2] + "b" #Format 2 letters + b 
                                    for atom_new in self.system_nopbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom:
                                        atom_new.residue_name = new_res_name             
                                    #######################################################################################################################
                                    crosslink_info.write("#Removed Branch Atoms. Renamed residue: " +  old_res_name + ".    New residue name: " + new_res_name + "\n")
                                    #######################################################################################################################                      
                                
                                #Delete the removable atom if specified for closing atom
                                elif t_sub_starting_closing_atoms[crosslinkNumber][1] != "null":    
                                    closing_idx = 0
                                    for atom_rmv in self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom:
                                        if atom_rmv.name == t_sub_starting_closing_atoms[crosslinkNumber][1]:
                                            atom_to_remove_close = atom_rmv #atom obect to be subsituted in the closing position (sub closing atom)
                                            removable_closing_idx = closing_idx
                                        closing_idx += 1                                   

                                    #Remove the removable atom from the closing residue
                                    del self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[removable_closing_idx]
                                    del self.system_nopbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[removable_closing_idx]

                                    #Rename the residue when atoms have been removed
                                    old_res_name = self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].name
                                    new_res_name = old_res_name[0:2] + atom_to_remove_close.name[-1]
                                    
                                    for atom_new in self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom:
                                        atom_new.residue_name = new_res_name
                                    for atom_new in self.system_nopbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom:
                                        atom_new.residue_name = new_res_name

                                    #######################################################################################################################
                                    crosslink_info.write("#Removed Atom: " + atom_to_remove_close.name + ". Renamed residue: " +  old_res_name +  ".    New residue name: " + new_res_name + "\n")
                                    #######################################################################################################################                                                                                            
                                
                                crosslink_info.write("bond    mol." + str(res_end_idx_tleap) + "." + atom_end.name + "    mol." + str(CL_end_idx_tleap) + "." + CL_atom_end.name + "\n")
                                crosslink_info.write("\n")

                                #Remove the molecule and the residue of starting and closing atoms of the [molecule,residue] so that they are not consider again for new crosslinks
                                dic_mol_res[starting_atom_name].remove([rnd_mol_ini, rnd_res_ini])
                                dic_mol_res[closing_atom_name].remove([new_mol_end, new_res_end])

                                #Process the next crosslinker molecule                
                                crosslinkNumber += 1        
                                self.system_nopbc.state[0] = len(t_molecules) + crosslinkNumber                    
                                self.system_pbc.state[0] = self.system_nopbc.state[0]
                                ##################################################
                            else:
                                self.bar.substract(len(self.system_pbc.t_molecules[len(t_molecules)+crosslinkNumber].t_residues))
                                del self.system_pbc.t_molecules[-1]
                                del self.system_nopbc.t_molecules[-1]
                                
                                if t_sub_starting_closing_atoms[crosslinkNumber][0] == "branch":
                                    for elem in t_removable_starting:
                                        self.coll.addAtomToDictionary(elem)
                                
                                elif t_sub_starting_closing_atoms[crosslinkNumber][0] != "null":
                                        self.coll.addAtomToDictionary(atom_to_remove_start)
                                
                        #########################################END BUILDING CROSSLINKS#######################################################        
                        
                        else:
                            t_crosslinkRetries[crosslinkNumber] = self.system_pbc.max_retries*tolerance + 1
                            t_errors.append('043') #No starting-closing atoms within the crosslink distance 


                    else: #No retries reminding
                        t_crosslinkRetries[crosslinkNumber] = 0
                        crosslinkNumber -= 1
                        
                        printRestart = True
                        t_errors.append('032') #Crosslinks couldn't be created
                        self.system_nopbc.state[0] = len(t_molecules)+crosslinkNumber 
                        self.system_pbc.state[0] = self.system_nopbc.state[0]
                        self.system_nopbc.state[2]=1
                        self.system_pbc.state[2]=1
                        v_exit       = True
                        crosslink_info.close()

            #Ending program
                crosslink_info.close()
            self.endSusi(t_errors,imTest,periodical_coordinates,language,box,minimumDistance, \
                         retries,maximum_deviation,pathPrepi,t_molecules,t_allResidues,t_pathFilePrepi,printRestart)

######################################################################################################################
######################################################################################################################
    def endSusi(self,t_errors,imTest,periodical_coordinates,language,box,minimumDistance,       \
                retries,maximum_deviation,pathPrepi,t_molecules,t_allResidues,t_pathFilePrepi,printRestart):
        """Check for error and print them. Close spacebar, save output files and print some information


        Input-----------------------
            t_errors (list): If some errors occured, this list contain all of them
            imTest (str): if equal to 'X' --> test mode activated, file.pdb will not be created
            periodical_coordinates (str): Use boundary condition / N: Don't use boundary condition
            language (str): Message language
            box (list): Array with 3 box dimensions (X,Y,Z) and unit
            minimumDistance (float): Size and unit of minimum distance lineage
            retries (int): Maximum number of retries to find a residue angle that don't collide with other atoms
            maximum_deviation (int): Value for maximum deviation
            pathPrepi (str): Path of PREPI files
            t_molecules (list): List of molecules obtained from file
            t_allResidues (list): List of all residues obtained from file
            t_pathFilePrepi (str): Path to find Prepi input files
            printRestart (bool): Whether printing a PDB restart file
        """

        try:
              # Terminate progress bar
            noErrors = True
            self.bar.terminate()
            if len(t_errors) == 0:
                t_errors.append('024')                #No error found
                t_errors.append('033')                #System written succesfully
                outputName1 = self.system_pbc.name+'.pdb'
                outputName2 = self.system_nopbc.name+'.pdb'
            else:
                noErrors = False
                self.printErrors(t_errors, language)  #Printing errors
                outputName2 = self.system_nopbc.name+'_restart.pdb'

            if imTest == '':
                if printRestart:
                    files.storeFile(outputName2, self.system_nopbc, language)
                elif noErrors: # Normal printing of system coordinates 
                    if periodical_coordinates == 'Y':
                        files.storeFile(outputName1, self.system_pbc, language)
                    else:
                        files.storeFile(outputName2, self.system_nopbc, language)
                else:
                    print("System empty.")

            else:
                print("System name:            " + self.system_pbc.name)
                print("X dimension size:       " + str(box[0]) + " " + box[1])
                print("Y dimension size:       " + str(box[2]) + " " + box[3])
                print("Z dimension size:       " + str(box[4]) + " " + box[5])
                print("Periodical coordinates: " + periodical_coordinates)
                print("Minimum distance:       " + str(minimumDistance[0]) + " " + minimumDistance[1])
                print("Retries:                " + str(retries))
                print("Maximum deviation:      " + str(maximum_deviation) + " %")
                print("Molecules number:       " + str(len(self.system_pbc.t_molecules)))
                print("Path PREPI files:       " + pathPrepi)
                print("Molecules:              " + str(t_molecules))
                print("Residue list:           " + str(t_allResidues))
                print("prepi files ok:         " + str(t_pathFilePrepi))

                self.printErrors(t_errors, language)
                sys.exit()

        except IOError as errorF:
            print(errorF.message)
            t_errors.append(errorF.message)
            self.printErrors(t_errors, language)
            sys.exit()

######################################################################################################################
######################################################################################################################
    def __addNewMolecule(self, it_molecule, iv_moleculeNumber, iv_moleculeID, to_filePREPI,t_residues=None,cyclic=False): 
        """
        Adds a new molecule in the system.
        
        Input-----------------------
            it_molecule (list): list of residues in the molecule.
            iv_moleculeNumber (int): index number of the molecule considering the whole system as a reference.
            iv_moleculeID (str): letter representing the name of the molecule in the system.
            to_filePREPI (list): list of all the prepi files 
        Output----------------------
        v_error (int): flag variable showing if the crosslinker was created
        """


        t_xyz              = [] #list of coordinates of every atom in a residue inculding dummies???
        v_error            = 0  #Value of error flag variable 
        t_last3m           = [] #list of coordinates of last 3 M-atoms including dummies if necesary
        t_aux_last3m       = [] 
        t_aux_xyz          = []
        t_history_last3m   = [0]
        t_history_xyz      = [0]
        t_residuRetries    = []  #Number of retries of residue number
        collision          = False
# 1-Create new molecule with molecule name from file. Here the molecule object is added with no processed residues
        v_molecule = molecule(iv_moleculeID, iv_moleculeID)
        self.system_pbc.addMolecule(v_molecule)
        v_molecule2 = molecule(iv_moleculeID, iv_moleculeID)
        self.system_nopbc.addMolecule(v_molecule2)
# ---------------------------------------------------------------------------------------------------------------------
        for fillResiduRetries in range(0,len(it_molecule)): t_residuRetries.append(0)
        vNumResidue = 0
        while vNumResidue < len(it_molecule):
            t_residuRetries[vNumResidue] += 1
            if self.verbosity >=2:                                       # Checking verbosity to avoid much time consuming writing sentence
                self.printLogs(' Residue '+ it_molecule[vNumResidue] + '(' + str(vNumResidue) + ')' + \
                        ' (retry:' + str(t_residuRetries[vNumResidue]) + ')' + \
                        '-'.join([str(i) for i in t_residuRetries]),2)
            if t_residuRetries[vNumResidue] < self.system_pbc.max_retries:   #Check if there are still retries for this residue
                t_aux_last3m = copy.deepcopy(t_history_last3m[vNumResidue]) #For 1st try & 1st residue is 0. For 1st residue is 0.
                t_aux_xyz    = copy.deepcopy(t_history_xyz) #1st try & 1st residue is [0]
                t_xyz, v_filePREPI, t_last3m = builder.__processResidue(self, vNumResidue, it_molecule[vNumResidue], \
                                                                      t_residuRetries[vNumResidue], iv_moleculeNumber, \
                                                                      to_filePREPI, t_aux_xyz, t_aux_last3m,t_residues,cyclic)
                #Molecules objects were created and added to the system before __processResidue. 
                #After processResidue, a list with the coordinates of all atoms of the residue is recovered, but no resiude has created yet
                
                collision = self.__addNewResidue(t_xyz, v_filePREPI, vNumResidue, iv_moleculeID, iv_moleculeNumber)
                
                #In __addNewResidue all the atoms objects were created. If no collisions are detected, the resiudes are created and placed in the system in the same function.
                # Then, if this happens the output of funtcion will be False. 
                # Otherwise, residues are neither created nor placed. And the retries will increment for the current residue

                if collision == False:                              #No collision, process next residue
                    vNumResidue +=1
                    self.system_pbc.state[1] = vNumResidue
                    self.system_nopbc.state[1] = self.system_pbc.state[1]
                    self.bar.add(1)                                  #Update progress bar counter
                    t_history_last3m.append(copy.deepcopy(t_last3m)) # [0, [[xm1,ym1,zm1],...,[xm3,ym3,zm3]],...]
                    t_history_xyz.append(copy.deepcopy(t_xyz))   # [0, [[x1,y1,z1],...,[xn,yn,zn]],...]

            else:
                if vNumResidue == 0:                                #No residue could be placed for this molecule. The molecule could not be built
                    return self.system_pbc, self.system_nopbc, 1
                else:
                    self.printLogs( '  Get previous residue', 2 )
                    t_residuRetries[vNumResidue]    = 0              #Reset retries for this residue
                    vNumResidue -=1                                  #Process again previous residue
                    self.bar.substract(1)                            #Update progress bar counter
                    #All information of the previous residue that was already placed in the system is deleted, so that it will be created again
                    del t_history_last3m[-1]
                    del t_history_xyz[-1]
                    self.coll.removeFromDictionary(self.system_pbc.t_molecules[iv_moleculeNumber].t_residues[-1].t_atom)
                    del self.system_pbc.t_molecules[iv_moleculeNumber].t_residues[-1]
                    del self.system_nopbc.t_molecules[iv_moleculeNumber].t_residues[-1]
# ---------------------------------------------------------------------------------------------------------------------
#v_error will be always 0 unless the program enters to the condition of residue num=0 and retries=maxretries        
        return v_error
######################################################################################################################
######################################################################################################################
    def __processResidue(self, iv_numResidue, iv_residueName, iv_numRetry, iv_moleculeNumber, it_filePREPI, \
                         it_xyz, it_last3m,t_residues=None,cyclic =False):
        """Get the coordinates of a residue in a molecule generated from the geometry of the prepi file.
        Input-----------------------
            iv_numResidue (int): number of the current residue within the molecule
            iv_residueName (str): name of the current residue
            iv_numRetry (int): number of retries for the current residue
            iv_moleculeNumber (int): index number of the molecule considering the whole system as a reference.
            it_filePREPI (list): list of all the prepi files
            it_xyz (list): list of the coordinates of all the residues placed in the system
            it_last3m (list): list of the coordinates of the last 3-M atoms of  the last residue placed for the crosslinker molecule (void if it is the first residue of the molecule)
        Output----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue 
            v_filePREPI (PREPI): name of the PREPI file for the current residue
            it_last3m (list) : list of the coordinates of the last 3-M atoms of the current residue
        """

        for x in range(0, len(it_filePREPI)):                     #Get correct filePREPI value from all filePREPI files
            if it_filePREPI[x].name == iv_residueName: 
                v_filePREPI = it_filePREPI[x]

        if iv_numRetry > 1 and iv_numResidue > 0:   #The second and next attempts to place the second and next residues. Modifies randomly th,ph from original PREPI info
            t_aux_par = builder.getRandomThetaAndPhi(v_filePREPI.t_par, self.system_pbc.maximum_deviation)
        else: #The first attempt to place the residue, from the second residue, uses r,th,ph from prepi File
            t_aux_par = v_filePREPI.t_par #list with r,angle,dihedral for every atom in PREPI including DUMMY
        if iv_numResidue == 0: #The first residue of the molecule always places the 1st of the 3 DUMMY on the origin
            it_xyz = [] #For the first residue is an empty list, so the residues are placed considering the dummy atoms as reference
            it_xyz    = int2car.int2car(v_filePREPI.t_con, t_aux_par, it_xyz) #it_xyz includes the coordiantes of all atoms of the resiude also 3 DUMMY. This list is initialized here always. The first dummy atom is placed on the origin.
            it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.

            if iv_moleculeNumber != 0:
            # First residue from second (and next) molecules must have random values for: Theta in 3 DUMM line and Theta and Phi in first M line after DUMM lines
                rX, rY, rZ, errorRandom = self.__getRandomPositionWithoutCollision()
                if errorRandom == True:
                    return it_xyz, v_filePREPI, it_last3m
                else:
                    for atom_in_res in range(0, len(it_xyz)): # From no DUMM lines (0,1,2) to all the atoms of the residue
                        it_xyz[atom_in_res][0] = it_xyz[atom_in_res][0] + rX
                        it_xyz[atom_in_res][1] = it_xyz[atom_in_res][1] + rY
                        it_xyz[atom_in_res][2] = it_xyz[atom_in_res][2] + rZ
                    it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
        else:  #iv_numResidue > 0
            if cyclic:
                idx_M0_p = self.__get1atomResidueName(t_residues[iv_numResidue-1],it_filePREPI)
                idx_M1_p = self.__getm1atomResidueName(t_residues[iv_numResidue-1],it_filePREPI)
                M0_p = it_xyz[iv_numResidue][idx_M0_p] #first M atom of the previous residue
                Mn_p = it_xyz[iv_numResidue][idx_M1_p]

            t_old_last3m = copy.deepcopy (it_last3m) #it_last3m always is initialized for the first residue of every molecule. The second residue will use info of M atoms of first residue if 2nd residue has less than 3 M-atoms
            it_xyz       = int2car.int2car(v_filePREPI.t_con, t_aux_par, it_last3m)#gets the coordinates using the modified th,ph if necessary. Uses last 3 M-atoms as reference no dummies. It places only no dummies atoms
            it_last3m    = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, t_old_last3m)#updates the last 3 M-atoms coordintaes, considering the info of the former residue. [[xM1,yM1,zM1],..,[xM3,yM3,zM3]]

            if cyclic:
                Vp = Mn_p - M0_p
                Vp /= np.linalg.norm(Vp)
                #extract the growth vector of the current residue
                idx_M0_n = self.__get1atomResidueName(iv_residueName,it_filePREPI)
                idx_M1_n = self.__getm1atomResidueName(iv_residueName,it_filePREPI)

                M0_n = it_xyz[idx_M0_n] #first M atom of the current residue
                Mn_n = it_xyz[idx_M1_n] #last M atom of the current residue

                Vn = Mn_n - M0_n 
                Vn /= np.linalg.norm(Vn)

                #establish a new plane with M0_p, Mn_p and  M0_n

                Vpn = M0_n - M0_p
                Vpn /= np.linalg.norm(Vpn)

                #Find a normal vector to the plane

                ax_N = np.cross(Vp,Vpn)
                ax_N /= np.linalg.norm(ax_N)

                #Find a rotation axis for Vn

                ax_K = np.cross(Vn,ax_N)
                ax_K /= np.linalg.norm(ax_K)

                #Find the angle of projection to the plane

                angle_plain = np.arccos(np.dot(Vn,ax_N)/(np.linalg.norm(Vn)*np.linalg.norm(ax_N)))
                angle_plain = np.pi/2 - angle_plain

                #Rotation of Vn to place it on the plane
                it_xyz = self.__Rot_Trans_Residue(it_xyz, ax_K, angle_plain,3, [M0_n[0],M0_n[1],M0_n[2]], None, False, 0)

                #Find the rotated point Mn_n 
                
                Mnr_n = it_xyz[idx_M1_n] #last M atom of the current residue after rotation to the plane
                Vnr = Mnr_n - M0_n 
                Vnr /= np.linalg.norm(Vnr)

                #Find the alignment angle of between the current and the previous residue

                angle_align = np.arccos(np.dot(Vp,Vnr))
                #Define the cyclic angle

                nres = len(t_residues)
                angle_cycle = 2*np.pi/(nres)
                angle_rot = angle_align - angle_cycle
                
                #Rotate the residue to form a cycle
                
                it_align_x = [list(vec) for vec in it_xyz]            
                it_xyz = self.__Rot_Trans_Residue(it_align_x, ax_N, angle_rot,3, [M0_n[0],M0_n[1],M0_n[2]], None, False, 0)
                #####Get the last 3M rotated atoms
                t_old_last3m = [it_xyz[0],it_xyz[1],it_xyz[2]]

                it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, t_old_last3m)#updates the last 3 M-atoms coordintaes, considering the info of the former residue.                                     
               
        return it_xyz, v_filePREPI, it_last3m
######################################################################################################################
######################################################################################################################
    def __addNewResidue(self, it_xyz, iv_filePREPI, iv_numResidue, iv_moleculeID, iv_moleculeNumber):
        """Check the collisions between the atoms of the new residue and all the existing atoms in the system. If there are no collisions the residue is added to the system.
        
        Input-----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue 
            iv_filePREPI (PREPI): name of the PREPI file for the current residue
            iv_numResidue (int): number of the current residue within the molecule
            iv_moleculeID (str): letter representing the name of the molecule in the system.
            iv_moleculeNumber (int): index number of the molecule considering the whole system as a reference.
        Output----------------------
            flag (bool): Boolean value. True if there were collsions in the system, false if there were not.
        """

# 1-Convert values obtained with function int2car to values in atom PDB format
        it_xyz2      = builder.convertxyz(it_xyz)
        t_aux_atoms  = builder.convertToAtoms(it_xyz2,iv_filePREPI.t_name, iv_filePREPI.name, iv_numResidue, iv_moleculeID, iv_filePREPI.t_topologic,iv_filePREPI.t_numbering,iv_filePREPI.t_con)
        t_aux_atoms2 = copy.deepcopy (t_aux_atoms)




# 2-Check if atoms are inside or outside box and calculate Periodic Boundary Condition if it's necessary
        for natoms in range(0, len(t_aux_atoms)):
            t_aux_atoms[natoms].x, t_aux_atoms[natoms].y, t_aux_atoms[natoms].z = \
                               builder.periodicBoundaryCondition(self.system_pbc.box, t_aux_atoms[natoms].x, \
                                                                                t_aux_atoms[natoms].y, \
                                                                                t_aux_atoms[natoms].z )
        if (iv_moleculeNumber)<0: #After the first molecule is placed it is necesary to study the posible growing paths
            t_aux_atoms = builder.__smartGrowth(self,t_aux_atoms,iv_filePREPI,self.system_pbc.maximum_deviation)
        
# 3-Check if any atom of new residue makes a collision with existing atoms
        if len(self.system_pbc.t_molecules) == 1 and len(self.system_pbc.t_molecules[0].t_residues) >= 1 or \
            ( len(self.system_pbc.t_molecules) > 1 ):
            for natoms in t_aux_atoms:    
                if natoms.name != 'DUMM':
                    if(self.coll.checkCollisionD(natoms)==True):
                        return True





# 4-Create residue with atoms and other variables
        v_residue  = residue(t_aux_atoms, iv_filePREPI.name, iv_numResidue, iv_moleculeID)
        v_residue2 = residue(t_aux_atoms2, iv_filePREPI.name, iv_numResidue, iv_moleculeID)
# 5-Delete DUMM lines from atom list (Don't want this lines in file that will be created later)
        v_residue.t_atom  = atom.deleteDUMMlines(v_residue.t_atom)
        for elem in v_residue.t_atom:
            self.coll.addAtomToDictionary(elem)
        v_residue2.t_atom = atom.deleteDUMMlines(v_residue2.t_atom)
# 6-Add residue to current system
        self.system_pbc.t_molecules[iv_moleculeNumber].addResidue(v_residue)
        self.system_nopbc.t_molecules[iv_moleculeNumber].addResidue(v_residue2)
        return False
######################################################################################################################
######################################################################################################################
    @staticmethod
    def convertxyz(t_xyz):
        """Converts value generated by int2car to values in format number of 8 positions with 3 decimals: %8.3f

        Input-----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue 
        Output----------------------
            t_auxxyz (list): list with the coordinates of all the atoms of the current residue in format number of 8 positions with 3 decimals:
        """
        t_newxyz = []
        for i in range(0,len(t_xyz)): #0<i < Number of atoms in t_xyz
            for j in range (0,len(t_xyz[i])):#j =0,1,2 -> x,y,z of atom i 
                t_auxxyz = builder.convertxyz2(t_xyz[i])
            t_newxyz.append(t_auxxyz)
        return t_newxyz
    @staticmethod
    def convertxyz2(t_xyz):
        t_auxxyz = []
        t_auxxyz.append( '%8.3f' % t_xyz[0])
        t_auxxyz.append( '%8.3f' % t_xyz[1])
        t_auxxyz.append( '%8.3f' % t_xyz[2])
        return t_auxxyz
######################################################################################################################
######################################################################################################################
    @staticmethod
    def convertToAtoms(t_xyz, t_atom_name, residueName, residueNumber, chainID, topology, numbering,connectivity):
        """Convert data obtained from conversion from class int2car to data in PDB format
        
        Input-----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue.
            t_atom_name (list): list with atom names.
            residueName (str): name of the residue.
            residueNumber (int): number index of the residue.
            chainID (str): letter representing the name of the molecule.
            topology (list): list with letters representing the toplogic type of an atom according to PREP format- 
            numbering (list): list with integers representing the indexing of atoms according to PREP format- ,
            connectivity (list): list of lists of integers representing the connectivity of atoms according to PREP format- ,
        
        Output----------------------
            t_atoms (list): list of atom objects.
        """

        t_atoms        = []
        v_tipus        = 'ATOM'
        v_serialNumber = 1
        v_occupancy    = 1.00
        v_temperature  = 0
        v_charg        = ''
        for i in range(0,len(t_atom_name)):
            v_xyz  = t_xyz[i]   #Position 0: X, Position 1: Y, Position 2: Z
#                                                                                 X         Y         Z
            v_atom = atom(v_tipus, v_serialNumber, t_atom_name[i], chainID, v_xyz[0], v_xyz[1], v_xyz[2], \
                          v_occupancy, v_temperature, t_atom_name[i][0], v_charg, residueName, residueNumber, topology[i], numbering[i], connectivity[i])
            t_atoms.append(v_atom)
      
            v_serialNumber +=1
        return t_atoms
######################################################################################################################
######################################################################################################################
    @staticmethod
    def printErrors(t_errors, language):
        """Print errors texts stored in a list
        
        Input-----------------------
        t_errors (list): list with error messages.
        language (str): Language for messages. English.
        """
        for n in range(0, len(t_errors)):
            if len(t_errors[n][0]) == 1:
                print(messages.get_text_message(t_errors[n], language))
            else:
                error = messages.get_text_message(t_errors[n][0], language)
                for numErrors in range(1, len(t_errors[n])):
                    error = str.replace(error, '&', t_errors[n][numErrors], 1)
                print(error)
######################################################################################################################
######################################################################################################################
    def printLogs(self,t_logs, verbLebel):
        """ Print log texts stored in a list
        Input-----------------------
            t_logs (list): List with logs.
            verbLevel (int): Level of verbosity allowed to be printed.

        """
        if verbLebel <= (self.verbosity):
            self.logFile.write(t_logs +'\n')
######################################################################################################################
######################################################################################################################
    @staticmethod
    def getLast3M(mLines, t_xyz, numResidue, t_old_last3m):
        """Search in a list last 3 values with M value and put in order.

        Input-----------------------
            mLines (int): Number value of lines with M value.
            t_xyz (list): list with the coordinates of all the atoms of the current residue. 
            numResidue (int): number of the current residue within the molecule.
            t_old_last3m (list): list of last 3M from previous residue.
        Output----------------------
        t_last3M (list): list of coordinates of the last 3 M-type atoms.
        """

        t_last3M = []
        v_mlines = len(mLines)  #Includes 3M Dumm lines
        if v_mlines   < 3: print ("Error in getLast3M")
        else:
            v_mlines -= 3 #Counts only the total number of M lines considering no DUMMY
# 4M = 1M = 2,3 DUMM lines + 1 M
# 5M = 2M =   3 DUMM line  + 2 M
            if v_mlines ==  1:
                if numResidue == 0: #First residue
#2nd DUMM line, 3th DUMM line, 1st M line not from DUMM lines
                    for v_mlines in range(1, 4): t_last3M.append(t_xyz[mLines[v_mlines]])
                else:
                    t_last3M.append(t_old_last3m[mLines[1]])  #2nd value from last3M from previous residue
                    t_last3M.append(t_old_last3m[mLines[2]])  #3th value from last3M from previous residue
                    t_last3M.append(t_xyz[mLines[3]])         #1st M line not from DUMM lines
            elif v_mlines ==  2:
                if numResidue == 0:
# 3th DUMM line, 1st M line not from DUMM lines, 2nd line not from DUMM lines
#                 t_last3M.append([0., 0., 0.])
                    for v_mlines in range(2, 5): t_last3M.append(t_xyz[mLines[v_mlines]])
                else:
                    t_last3M.append(t_old_last3m[mLines[2]])  #3th value from last3M from previous residue
# 1st M line not from DUMM lines, 2nd M line not from DUMM lines
                    for v_mlines in range(3, 5): t_last3M.append(t_xyz[mLines[v_mlines]])
            else:
                for v_mlines in range(v_mlines,len(mLines)):
                    t_last3M.append(t_xyz[mLines[v_mlines]])
        return t_last3M
######################################################################################################################
######################################################################################################################
    @staticmethod
    def periodicBoundaryCondition(box, posX, posY, posZ):
        """Checks if any dimension of atoms position it's outside box, if this happens calculate Periodic Boundary Condition
        box (list): Box dimensions. box[0]: Dimension X. box[2]: Dimension Y. box[4]: Dimension Z.

        Input-----------------------
        posX (float): X-axis coordinate.
        posY (float): Y-axis coordinate.
        posZ (float): Z-axis coordinate.
        Output----------------------
        posX (float): X-axis coordinate with PBC repositioning.
        posY (float): Y-axis coordinate with PBC repositioning.
        posZ (float): Z-axis coordinate with PBC repositioning.
        """

        posX = posX - box[0]*np.rint(posX/box[0])
        posY = posY - box[2]*np.rint(posY/box[2])
        posZ = posZ - box[4]*np.rint(posZ/box[4])
        return posX, posY, posZ
######################################################################################################################
######################################################################################################################
    def __checkCollisionInSystem2(self, posX, posY, posZ):
        """Check if XYZ have a collision with any atom in the system.
        
        Input-----------------------
        posX (float): X-axis coordinate.
        posY (float): Y-axis coordinate.
        posZ (float): Z-axis coordinate.
        Output----------------------
        flag (bool): True if a collision is detected.
        """

        p1 = np.array([posX, posY, posZ])
        for nMolecules in range(0,len(self.system_pbc.t_molecules)):
            for nResidues in range (0, len(self.system_pbc.t_molecules[nMolecules].t_residues)):
                for nAtoms in range (0,len(self.system_pbc.t_molecules[nMolecules].t_residues[nResidues].t_atom)):
                    p2       = np.array([self.system_pbc.t_molecules[nMolecules].t_residues[nResidues].t_atom[nAtoms].x, \
                                         self.system_pbc.t_molecules[nMolecules].t_residues[nResidues].t_atom[nAtoms].y, \
                                         self.system_pbc.t_molecules[nMolecules].t_residues[nResidues].t_atom[nAtoms].z])
                    distance = np.linalg.norm(p1-p2)
                    if distance < self.system_pbc.minimum_distance[0]:
                        return True
        return False
#######################################################################################################################
######################################################################################################################
    def __getRandomPositionWithoutCollision(self):
        """Try to get random position N times. If no position can't be found, then returns an error.

        Output----------------------
        posX (float): X-axis coordinate.
        posY (float): Y-axis coordinate.
        posZ (float): Z-axis coordinate.
        flag (bool): True if a collision is detected.
        """
        #RANDOMSEED
        #random.seed(1)


        for num in range(1, int(self.system_pbc.max_retries)):
            x = random.uniform(0, self.system_pbc.box[0])
            y = random.uniform(0, self.system_pbc.box[2])
            z = random.uniform(0, self.system_pbc.box[4])
            if self.__checkCollisionInSystem2(x, y, z) ==  False:
                return x, y, z, False
        return x, y, z, True
######################################################################################################################
######################################################################################################################
    @staticmethod
    def getRandomThetaAndPhi(t_par, maximum_deviation):
        """Read values from list t_par, and change value of Theta in 3 DUMM line and Theta and Phi in first M line after DUMM lines.
        THETA  it's the bond angle between atom NB(I), NA(I) and I.
        PHI(I) it's the dihedral angle between NC(I), NB(I), NA(I) and I.

        Input-----------------------
        t_par (list): t_par lines from PREPI file containing the information of internal coordinates.
        maximum_deviation (float): Value for maximum deviation to deformate internal coordinates.
        Output----------------------
        t_new_par (list): t_par lines from PREPI file containing the information of internal coordinates deformated.
        """    

        t_new_par = copy.deepcopy(t_par)
        if len(t_new_par) >= 1:
# We induce distortion in all residue parameters
            #RANDOMSEED
            #random.seed(1)
            for num in range(0,len(t_new_par)):
                for ele in range(1,3): #Changes only the angle and th dihedral #original 0,2 changes only distance and angle
                    randomValue = ( 10000 - float(random.randrange(int(-maximum_deviation*100), int(maximum_deviation*100)))) / 10000
                    t_new_par[num][ele]           = t_new_par[num][ele] * randomValue      # r
        return t_new_par
######################################################################################################################
######################################################################################################################
    def getSystem(self):
        """
        Returns the created system with pbc
        """
        return self.system_pbc
######################################################################################################################
######################################################################################################################
    def RotMat(angle,director):
        """
        Rotation matrixes.

        Input-----------------------
            angle (float): angle in radians
            director (list): untary director vector
        
        Output----------------------
            matR (list): Rotation matrix
        """
        theta = angle
        u = director

        #theta = np.pi*theta/180.0 only if angle is in degrees
        ux = u[0]; uy = u[1]; uz = u[2]
        s = np.sin(theta); c = np.cos(theta) ; mc = 1 - c

        matR = [[c+ux**2*mc, ux*uy*mc-uz*s, ux*uz*mc+uy*s],
               [uy*ux*mc+uz*s, c+uy**2*mc, uy*uz*mc-ux*s],
               [uz*ux*mc-uy*s, uz*uy*mc+ux*s, c+uz**2*mc]]
        
        return matR
######################################################################################################################
######################################################################################################################
    def __rot_arb(self,point,axis,angle):
        """
        Rotation around an arbitrary axis basis on fractal trees growth.

        Input-----------------------
            point (list): coordinates of the point to be rotated.
            axis (list): vector representing the axis around which the point will be rotated.
            angle (float): rotation angle.
        Output----------------------
            rotate_cw (list): clockwise rotated point
            rotate_ccw (list):vcounter clockwise rotated point
        """
        #P' = (n·P)*n + cos(theta)*[P - (n·P)*n] +- sin(theta)*(n X P)

        rotate_cw = np.dot(axis,point)*axis + np.cos(angle)*(point-np.dot(axis,point)*axis) + np.sin(angle)*(np.cross(axis,point))
        rotate_ccw = np.dot(axis,point)*axis + np.cos(angle)*(point-np.dot(axis,point)*axis) - np.sin(angle)*(np.cross(axis,point))
        
        return rotate_ccw
######################################################################################################################
######################################################################################################################
    def __smartGrowth(self,t_aux_atoms,iv_filePREPI,maximum_deviation):
        """
        Strategy using our 5% deviation every time we add a residue to build up our residue to the places in the space with lower density that is equal to less collision probability  

        Input-----------------------
            t_aux_atoms (list): The list with atoms in the residue-
            iv_filePREPI (filePREPI): The prepi file our residue: list of float
        Output----------------------
            t_aux_atoms (list): The list with the atoms modified
        
        """       
        firstMatom = iv_filePREPI.t_mLines[3]
        lastMatom  = iv_filePREPI.t_mLines[-1]

        vector=[t_aux_atoms[lastMatom].x-t_aux_atoms[firstMatom].x,t_aux_atoms[lastMatom].y-t_aux_atoms[firstMatom].y, t_aux_atoms[lastMatom].z-t_aux_atoms[firstMatom].z] #We generate the vector of the residue's growth
        vector1=[math.sqrt(vector[0]**2+vector[1]**2+vector[2]**2),math.atan2(vector[1],vector[0]),math.atan2(math.sqrt(vector[0]**2+vector[1]**2),vector[2])] #we transforme it into espheric coordinates
        vu=[0,0,0] #Initializing the residue's growth unitary vector
        vu[0]=vector[0]/vector1[0] #X coordinate of unitary vector
        vu[1]=vector[1]/vector1[0] #Y coordinate of unitary vector
        vu[2]=vector[2]/vector1[0] #Z coordinate of unitary vector
        vector1[2]=vector1[2]-math.pi*(maximum_deviation/100) # We apply the 5% deviation that our residue can have in one axis
        vector1=[vector1[0]*math.cos(vector1[1])*math.sin(vector1[2]),vector1[0]*math.sin(vector1[1])*math.sin(vector1[2]),vector1[0]*math.cos(vector[2])] #We recover against the X,Y,Z coordinates of the vector modified
        rotation_angle = np.pi/4
        Mrot = builder.RotMat(rotation_angle,vu)
        vector=[[t_aux_atoms[lastMatom].x-t_aux_atoms[firstMatom].x,t_aux_atoms[lastMatom].y-t_aux_atoms[firstMatom].y, t_aux_atoms[lastMatom].z-t_aux_atoms[firstMatom].z]] #Now we create the same vector as beginning but now in a list 
        vector.append(vector1) #we append the vector with deviation 
        for x in range(7): # This loop rotates vector1 every time pi/4 more than the last one to create all the possible positions
            vector2 = copy.deepcopy(vector[-1])           
            for y in range(3):
                vector2[y]=vector2[0]*Mrot[y][0]+vector2[1]*Mrot[y][1]+vector2[2]*Mrot[y][2]
            vector.append(vector2) # we append all the vectors in the list  
        density=100000000 #Initialization of the density with a high number
        savekey=[]
        for x in range(9): #This loop generates all the positions and calculates the density giving more value at shorter distances
            at=atom(0,0,0,0,t_aux_atoms[0].x+vector[x][0],t_aux_atoms[0].y+vector[x][1],t_aux_atoms[0].z+vector[x][2],0,0,0,0,0,0,0,0,0)
            at.x, at.y, at.z=builder.periodicBoundaryCondition(self.system_pbc.box,at.x,at.y,at.z)
            density1=3*self.coll.getLengthOfList(at)
            for y  in range(9):
                at2=atom(0,0,0,0,t_aux_atoms[0].x+vector[x][0]+vector[y][0],t_aux_atoms[0].y+vector[x][1]+vector[y][1],t_aux_atoms[0].z+vector[x][2]+vector[y][2],0,0,0,0,0,0,0,0,0)
                at2.x, at2.y, at2.z=builder.periodicBoundaryCondition(self.system_pbc.box,at2.x,at2.y,at2.z)
                density2=2*self.coll.getLengthOfList(at2)+density1
                for z in range(9):
                    at3=atom(0,0,0,0,t_aux_atoms[0].x+vector[x][0]+vector[y][0]+vector[z][0],t_aux_atoms[0].y+vector[x][1]+vector[y][1]+vector[z][1],t_aux_atoms[0].z+vector[x][2]+vector[y][2]+vector[z][2],0,0,0,0,0,0,0,0,0)
                    at3.x, at3.y, at3.z=builder.periodicBoundaryCondition(self.system_pbc.box,at3.x,at3.y,at3.z)
                    density3=self.coll.getLengthOfList(at3)+density2
                    if density3 < density: #At the end we take the lowest density and the first movement vector
                        savekey=[]
                        density=density3*1
                        savekey.append(x)
                    if density3 == density:
                        savekey.append(x)
        vectorkey=random.choice(savekey)
        vectr=[t_aux_atoms[lastMatom].x-t_aux_atoms[firstMatom].x,t_aux_atoms[lastMatom].y-t_aux_atoms[firstMatom].y, t_aux_atoms[lastMatom].z-t_aux_atoms[firstMatom].z] #We generate the vector of the residue's growth
        vector1=[math.sqrt(vectr[0]**2+vectr[1]**2+vectr[2]**2),math.atan2(vectr[1],vectr[0]),math.atan2(math.sqrt(vectr[0]**2+vectr[1]**2),vectr[2])] #Transformation in spherical coordinates
        vector2=vector[vectorkey]*1 #The best movement vector choose with the previous loop
        vector2=[math.sqrt(vector2[0]**2+vector2[1]**2+vector2[2]**2),math.atan2(vector2[1],vector2[0]),math.atan2(math.sqrt(vector2[0]**2+vector2[1]**2),vector2[2])] #Transformation in spherical coordinates
        diftheta=vector2[1]-vector1[1] #Difference of the first angle
        difpsi=vector2[2]-vector1[2] #Difference of the second angle
        vector=[t_aux_atoms[1].x-t_aux_atoms[0].x,t_aux_atoms[1].y-t_aux_atoms[0].y, t_aux_atoms[1].z-t_aux_atoms[0].z] # Now we generate the vector between the 2 first atoms
        vector1=[math.sqrt(vector[0]**2+vector[1]**2+vector[2]**2),math.atan2(vector[1],vector[0]),math.atan2(math.sqrt(vector[0]**2+vector[1]**2),vector[2])] #We transformate it in spherical coordinates
        vector1[1]=vector1[1]+diftheta # We apply the difference obtained above for the first angle
        vector1[2]=vector1[2]+difpsi #We apply the difference obtained above for the first angle 
        vector1=[vector1[0]*math.cos(vector1[1])*math.sin(vector1[2]),vector1[0]*math.sin(vector1[1])*math.sin(vector1[2]),vector1[0]*math.cos(vector1[2])] #We transform it into X,Y,Z coordinates
        difvector=[] #The list where we will save the difference between the atoms 
        llargada=len(t_aux_atoms)-1 #Parameter for the length of next loop
        for x in range (1,llargada): #Bucle for generating and save the difference between continuous atoms, it starts from 1 because between the first 2 atoms we apply the vector with the angle differences
            vec=[t_aux_atoms[x+1].x-t_aux_atoms[x].x,t_aux_atoms[x+1].y-t_aux_atoms[x].y, t_aux_atoms[x+1].z-t_aux_atoms[x].z]
            difvector.append(vec)
        t_aux_atoms[1].x=t_aux_atoms[0].x+vector1[0] #Application of the difference between the first 2 atoms coordinate X
        t_aux_atoms[1].y=t_aux_atoms[0].y+vector1[1] #Application of the difference between the first 2 atoms coordinate Y
        t_aux_atoms[1].z=t_aux_atoms[0].z+vector1[2] #Application of the difference between the first 2 atoms coordinate Z
        for x in range (1,llargada): #Bucle to create the new position between atoms once the 2nd atom is displaced conserving the differences
            v1=difvector[x-1]
            t_aux_atoms[x+1].x=t_aux_atoms[x].x+v1[0]
            t_aux_atoms[x+1].y=t_aux_atoms[x].y+v1[1]
            t_aux_atoms[x+1].z=t_aux_atoms[x].z+v1[2]
            t_aux_atoms[x+1].x,t_aux_atoms[x+1].y,t_aux_atoms[x+1].z=builder.periodicBoundaryCondition(self.system_pbc.box,t_aux_atoms[x+1].x,t_aux_atoms[x+1].y,t_aux_atoms[x+1].z) #We apply boundary conditions 
        return t_aux_atoms
######################################################################################################################
######################################################################################################################
    def __get1atomResidue(self,residue,to_filePREPI):
        """
        Get the first M-type atom of a residue
        
        Input-----------------------
            residue (Residue): Residue type object
            to_filePREPI (list): list of filePREPI
        Output----------------------
            firstatom (int): index of the first M-type atom of the residue

        """
        for x in range(len(to_filePREPI)):
            if residue.name == to_filePREPI[x].name: v_filePREPI=to_filePREPI[x]
        firstatom = v_filePREPI.t_mLines[3] 
        return firstatom
######################################################################################################################
######################################################################################################################
    def __get1atomResidueName(self,residue,to_filePREPI):
        """
        Get the first M-type atom of a residue
        
        Input-----------------------
            residue (str): name of the residue
            to_filePREPI (list): list of filePREPI
        Output----------------------
            firstatom (int): index of the first M-type atom of the residue

        """
        for x in range(len(to_filePREPI)):
            if residue == to_filePREPI[x].name: v_filePREPI=to_filePREPI[x]
        firstatom = v_filePREPI.t_mLines[3] 
        return firstatom
######################################################################################################################
######################################################################################################################
    def __getm1atomResidue(self,residue,to_filePREPI):
        """
        Get the last M-type atom of a residue
        
        Input-----------------------
            residue (Residue): Residue type object
            to_filePREPI (list): list of filePREPI
        Output----------------------
            m1atom (int): index of the last M-type atom of the residue
            
        """
        for x in range(len(to_filePREPI)):
            if residue.name == to_filePREPI[x].name: v_filePREPI=to_filePREPI[x]
        m1atom = v_filePREPI.t_mLines[-1] 
        return m1atom
######################################################################################################################
######################################################################################################################
    def __getm1atomResidueName(self,residue,to_filePREPI):
        """
        Get the last M-type atom of a residue
        
        Input-----------------------
            residue (str): name of the residue
            to_filePREPI (list): list of filePREPI
        Output----------------------
            m1atom (int): index of the last M-type atom of the residue
            
        """
        for x in range(len(to_filePREPI)):
            if residue == to_filePREPI[x].name: v_filePREPI=to_filePREPI[x]
        m1atom = v_filePREPI.t_mLines[-1] 
        return m1atom

######################################################################################################################
######################################################################################################################
    def __addNewCrossLink(self, it_molecule, iv_moleculeNumber, iv_moleculeID, to_filePREPI,vector_growth,vector_end_ini,atom_ini,atom_end_idx, atom_end_name,t_res_cl, iv_totalResidues,it_mol_res,iv_mol_ini,iv_distortion_angle, iv_bond_length_CL_start, iv_bond_length_CL_close,sub_closing_name):
        """
        Adds a new crosslinker molecule in the system.
        
        Input-----------------------
            it_molecule (list): list of residues in the crosslinker molecule.
            iv_moleculeNumber (int): index number of the crosslinker molecule considering the whole system as a reference.
            iv_moleculeID (str): letter representing the name of the crosslinker molecule in the system.
            to_filePREPI (list): list of all the prepi files
            vector_growth (list): vector of growth direction of the crosslinker molecule
            vector_end_ini (list): vector direction between the initial and final atoms to be substituted
            atom_ini (atom): starting atom object. 
            atom_end_idx (list): list with the indexes of the molecule, the residue and the closing atom.
            atom_end_name (str): name of the closing atom
            t_res_cl (list): list with residue objects of the crosslinker molecule. The first one has dummy atoms.
            iv_totalResidues (int): total number of residues of the molecule which the current residue belongs.
            it_mol_res (list): list containing the indexes of the [molecules,resiudes] in the system where closing-type atoms are found
            iv_mol_ini (int): index number of the molecule where the starting atom is located
            iv_distortion_angle (float): value of the maximum distortion angle a crosslink residue can rotate to avoid a colission and point to a new closing atom
            iv_bond_length_CL_start (float): Bond distance between the first atom of the crosslink and the starting atom
            iv_bond_length_CL_close (float): Bond distance between the last atom of the crosslink ande the closing atom
            sub_closing_name: name of the atom to be substituted to close the crosslink.     
        Output----------------------
        v_error (int): flag variable showing if the crosslinker was created
        """
        t_xyz              = [] #list of coordinates of every atom in a residue inculding dummies
        v_error            = 0  #Value of error flag variable 
        t_last3m           = [] #list of coordinates of last 3 M-atoms including dummies if necesary
        t_aux_last3m       = [] 
        t_aux_xyz          = []
        t_history_last3m   = [0]
        t_history_xyz      = [0]
        t_residuRetries    = [] #Number of retries of residue number
        collision          = False
        #iv_moleculeNumber is the current number counting both molecules and the current crosslink
# 1-Create new molecule with molecule name from file. Here the molecule object is added with no processed residues
        v_molecule = molecule(iv_moleculeID, iv_moleculeID)
        self.system_pbc.addMolecule(v_molecule)
        v_molecule2 = molecule(iv_moleculeID, iv_moleculeID)
        self.system_nopbc.addMolecule(v_molecule2)
# ---------------------------------------------------------------------------------------------------------------------
        for fillResiduRetries in range(0,len(it_molecule)): t_residuRetries.append(0)
        vNumResidue = 0
        while vNumResidue < len(it_molecule):

            t_residuRetries[vNumResidue] += 1

            if t_residuRetries[vNumResidue] < self.system_pbc.max_retries:   #Check if there are still retries for this residue
                t_aux_last3m = copy.deepcopy(t_history_last3m[vNumResidue])
                t_aux_xyz    = copy.deepcopy(t_history_xyz[vNumResidue])
                
                t_xyz, v_filePREPI, t_last3m, new_mol_end, new_res_end, new_HM_atom_end, change_flag = builder.__processResidueCrossLink(self, vNumResidue, it_molecule[vNumResidue],t_residuRetries[vNumResidue],\
                                                                      to_filePREPI, t_aux_xyz, t_aux_last3m,vector_growth,vector_end_ini, \
                                                                        atom_ini,atom_end_idx, atom_end_name,t_res_cl,it_mol_res,iv_mol_ini, iv_distortion_angle, iv_bond_length_CL_start, iv_bond_length_CL_close )
                
                #The molecule object of the new crosslink has been created
                #After __processResidueCrossLink the first residue has been placed in the ini atom postion, and the following residues (and atoms) have been created and placed.
                #The crosslink molecule object does not have any residues yet.
 
                if new_mol_end != None and new_res_end != None and new_HM_atom_end != None:
                    ####Removing from collision dictionary the heavy atoms when its specified
                    if sub_closing_name != "null" and sub_closing_name != "branch":    
                        for atom_rmv in self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom:
                            if atom_rmv.name == sub_closing_name:                               
                                atom_to_remove_close = atom_rmv #atom obect to be subsituted in the closing position (sub closing atom)
                                if self.coll.checkAtomsinList(atom_to_remove_close,atom_to_remove_close):                                   
                                    self.coll.removeFromDictionary([atom_to_remove_close])

                    elif sub_closing_name == "branch":
                        #Atoms list of the closing atom residue. Does not contain DUMMY atoms
                        t_atoms_closing = self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom
                        #Get indexes of atoms connected to the starting atom in the residue
                        t_connected_closing, in_closing_atom = self.__find_connected_atoms(t_atoms_closing,atom_end_name)
                        #Get a list of non-M heavy atoms connected to the starting atom
                        t_removable_closing = []
                        for i in t_connected_closing:
                            if t_atoms_closing[i].topology_type != 'M':
                                if 'H' not in t_atoms_closing[i].name:
                                    t_removable_closing.append(t_atoms_closing[i])
                                elif t_atoms_closing[i].connectivity[0] != in_closing_atom:
                                    t_removable_closing.append(t_atoms_closing[i])

                        for elem_in_rem in t_removable_closing:
                            if self.coll.checkAtomsinList(elem_in_rem,elem_in_rem):
                                self.coll.removeFromDictionary([elem_in_rem])

                    collision = self.__addNewResidueCrossLink(t_xyz, v_filePREPI, vNumResidue, iv_moleculeID, iv_moleculeNumber)

                    #In __addNewResidueCL all the atoms objects were created. If no collisions are detected, the resiudes are created and placed in the system in the same function.
                    # If no collisions are detected the output of the funtcion will be False. 
                    # Otherwise, residues are neither created nor placed. And the retries will increment for the current residue
                    if collision == False :                              #No collision, process next residue
                        vNumResidue +=1
                        self.system_pbc.state[1] = vNumResidue
                        self.system_nopbc.state[1] = self.system_pbc.state[1]
                        self.bar.add(1)                                  #Update progress bar counter
                        t_history_last3m.append(copy.deepcopy(t_last3m)) # [0, [[xm1,ym1,zm1],...,[xm3,ym3,zm3]],...]
                        t_history_xyz.append(copy.deepcopy(t_xyz))   # [0, [[x1,y1,z1],...,[xn,yn,zn]],...]
                        atom_end_idx[0] = new_mol_end
                        atom_end_idx[1] = new_res_end
                        atom_end_idx[2] = new_HM_atom_end
                                               
                    else:
                        ####Add to the collision dictionary the previously removed heavy atoms after a failled of placing the crosslink when heavy its specified
                        if sub_closing_name == "branch":
                            for elem in t_removable_closing:
                                if not(self.coll.checkAtomsinList(elem,elem)):
                                    self.coll.addAtomToDictionary(elem)
                                
                        elif sub_closing_name != "null": 
                            if not(self.coll.checkAtomsinList(atom_to_remove_close,atom_to_remove_close)):
                                self.coll.addAtomToDictionary(atom_to_remove_close)
                        
            else:
                if vNumResidue == 0:                #No residue could be placed for this molecule. The molecule could not be built
                    return 1, new_mol_end, new_res_end, new_HM_atom_end
                else:
                    self.printLogs( '  Get previous residue', 2 )
                    t_residuRetries[vNumResidue]    = 0              #Reset retries for this residue
                    vNumResidue -= 1                                  #Process again previous residue
                    self.bar.substract(1)                            #Update progress bar counter
                    #All information of the previous residue that was already placed in the system is deleted, so that it will be created again
                    del t_history_last3m[-1]
                    del t_history_xyz[-1]
                    self.coll.removeFromDictionary(self.system_pbc.t_molecules[iv_moleculeNumber].t_residues[-1].t_atom)
                    del self.system_pbc.t_molecules[iv_moleculeNumber].t_residues[-1]
                    del self.system_nopbc.t_molecules[iv_moleculeNumber].t_residues[-1]
# ---------------------------------------------------------------------------------------------------------------------
        return v_error, new_mol_end, new_res_end, new_HM_atom_end
######################################################################################################################
######################################################################################################################
    def __processResidueCrossLink(self, iv_numResidue, iv_residueName, iv_numRetry, it_filePREPI, \
                         it_xyz, it_last3m,vector_growth,vector_end_ini,atom_ini, atom_end_idx, atom_end_name,t_res_cl,t_mol_res,v_mol_ini, v_distortion_angle, v_bond_length_CL_start, v_bond_length_CL_close):
        """Get the coordinates of a residue in a molecule generated from the geometry of the prepi file.
        
        Input-----------------------
            iv_numResidue (int): number of the current residue within the crosslinker molecule
            iv_residueName (str): name of the current residue
            iv_numRetry (int): number of retries for the current residue 
            it_filePREPI (list): list of all the prepi files
            it_xyz (list): list of the coordinates of all the residues placed in the system
            it_last3m (list): list of the coordinates of the last 3-M atoms of  the last residue placed for the crosslinker molecule (void if it is the first residue of the molecule)
            vector_growth (list): vector of growth direction of the crosslinker molecule
            vector_end_ini (list): vector direction between the starting and closing atoms
            atom_ini (atom): starting atom. 
            atom_end_idx (list): list with the indexes of the molecule, the residue and the closing atom.
            atom_end_name (str): closing atom name.
            t_res_cl (list): list with residue objects of the crosslinker molecule. The first one has dummy atoms.
            t_mol_res (list): list containing the indexes of the [molecules,resiudes] in the system where closing-type atoms are found
            v_mol_ini (int): index number of the molecule where the starting atom is located
            v_distortion_angle (float): value of the maximum distortion angle a crosslink residue can rotate to avoid a colission and point to a new closing atom
            v_bond_length_CL_start (float): Bond distance between the first atom of the crosslink and the starting atom
            v_bond_length_CL_close (float): Bond distance between the last atom of the crosslink and the closing atom


        Output----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue.
            v_filePREPI (PREPI): name of the PREPI file for the current residue.
            it_last3m (list) : list of the coordinates of the last 3-M atoms of the current residue.
            v_error (int): flag variable showing if the residues were placed. 
            change_flag (bool): whether a new closing atom has been found.

        """
        for x in range(0, len(it_filePREPI)):                     #Get correct filePREPI value from all filePREPI files
            if it_filePREPI[x].name == iv_residueName: v_filePREPI = it_filePREPI[x]

        if iv_numRetry > 1: #The second and next attempts to place the next residues. Modifies randomly th,ph from original PREPI info
            t_aux_par = builder.getRandomThetaAndPhi(v_filePREPI.t_par, self.system_pbc.maximum_deviation)
        else: #The first attempt to place the residue, from the second residue, uses r,th,ph from prepi File
            t_aux_par = v_filePREPI.t_par #list with r,angle,dihedral for every atom in PREPI including DUMMY

        if iv_numResidue == 0: #The first residue of the molecule will be placed on the end-ini vector pointing the closing atom at a distince v_bond_length_CL_start from atom_ini      
            angle = np.arccos(np.dot(vector_end_ini,vector_growth)) #Angle between crosslink growth direction and end-ini substitution atoms direction
            axis = np.cross(vector_end_ini,vector_growth)# orthonormal crosslink growth direction and end-ini substitution atoms direction
            
            axis = axis/np.linalg.norm(axis)# unitary orthonormal crosslink growth direction and end-ini substitution atoms direction
            it_res_cl = copy.deepcopy(t_res_cl)
            it_res0_cl = [] #list that will containt the atoms objects of the first residue of crosslink
            
            for atom_cl in it_res_cl[0].t_atom:
                x_atom = atom_cl.x
                y_atom = atom_cl.y  
                z_atom = atom_cl.z  
                it_res0_cl.append([x_atom,y_atom,z_atom])
            
            #Deform the 1st residue if collisions were detected
            if iv_numRetry > 1:
                it_res0_cl = int2car.int2car(v_filePREPI.t_con, t_aux_par, [np.array(it_res0_cl[0]),np.array(it_res0_cl[1]),np.array(it_res0_cl[2])])

            it_xyz = self.__Rot_Trans_Residue(it_res0_cl, axis, angle,3, [atom_ini.x, atom_ini.y, atom_ini.z], vector_end_ini, True, v_bond_length_CL_start )#3 is the first no dummy M atom
            it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
                     

            new_mol_end = atom_end_idx[0]
            new_res_end = atom_end_idx[1]
            new_HM_atom_end = atom_end_idx[2]

            change_flag = False

        else:  #iv_numResidue > 0
            t_old_last3m = copy.deepcopy (it_last3m) #it_last3m always is initialized for the first residue of every molecule. The second residue will use info of M atoms of first residue if 2nd residue has less than 3 M-atoms
            it_xyz       = int2car.int2car(v_filePREPI.t_con, t_aux_par, it_last3m)#Gets the coordinates using the modified th,ph if necessary. Uses last 3 M-atoms as reference. 
                                                                                   #It places only no dummy atoms, the first 3 in the list are the last 3 M-atoms of the previous residue
                                                                                   #These first 3 M-atoms will be labeled as dummies by the convertToAtoms method
            it_last3m    = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, t_old_last3m)#updates the last 3 M-atoms coordintaes, considering the info of the former residue.
            new_mol_end = atom_end_idx[0]
            new_res_end = atom_end_idx[1]
            new_HM_atom_end = atom_end_idx[2]
            change_flag = False

            #if iv_numRetry > 1:
            if iv_numResidue < len(t_res_cl)-1 and iv_numRetry > 1:

                #Measure the size of the crosslinker starting from the current residue but growing from the next one
                dist_rm, atomM_last, atomM_ini, vector_growth_rm, _,_,_ = self.__ReMeasureMolecule( t_res_cl, iv_numResidue, v_filePREPI, it_filePREPI, it_last3m, it_xyz,True,False,False,False)
                    
                ###########Searching a new closing atom #############
                new_mol_end, new_res_end, new_HM_atom_end =  self.__get_atom_end_idx( t_mol_res, dist_rm, v_bond_length_CL_start,v_bond_length_CL_close, atomM_ini, v_mol_ini, atom_end_name, True)
                if new_mol_end != None:
                    #new_atom_end = self.system_nopbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[new_HM_atom_end]
                    new_atom_end = self.system_pbc.t_molecules[new_mol_end].t_residues[new_res_end].t_atom[new_HM_atom_end]
                    vector_end_ini_rm = np.array([new_atom_end.x, new_atom_end.y, new_atom_end.z]) - np.array([atomM_ini.x, atomM_ini.y, atomM_ini.z])
                    vector_end_ini_rm = vector_end_ini_rm/np.linalg.norm(vector_end_ini_rm)
                   
                    angle_rm = np.arccos(np.dot(vector_end_ini_rm,vector_growth_rm)) #Angle between crosslink growth direction and end-ini substitution atoms direction
                    
                    if angle_rm*180/np.pi > v_distortion_angle:
                        new_mol_end = None
                        new_res_end = None
                        new_HM_atom_end = None
                        change_flag = False
                    else:
                        axis_rm = np.cross(vector_end_ini_rm,vector_growth_rm)# orthonormal crosslink growth direction and end-ini substitution atoms direction
                        axis_rm = axis_rm/np.linalg.norm(axis_rm)# unitary orthonormal crosslink growth direction and end-ini substitution atoms direction
                        ######
                        it_xyz = self.__Rot_Trans_Residue(it_xyz, axis_rm, angle_rm,3, [atomM_ini.x, atomM_ini.y, atomM_ini.z], None, False, 0)              
                        #####Get the last 3M rotated atoms
                        t_old_last3m = [it_xyz[0],it_xyz[1],it_xyz[2]]
                        it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, t_old_last3m)#updates the last 3 M-atoms coordintaes, considering the info of the former residue.                                     
                        change_flag = True
                    
        return it_xyz, v_filePREPI, it_last3m, new_mol_end, new_res_end, new_HM_atom_end, change_flag
######################################################################################################################
######################################################################################################################
    def __addNewResidueCrossLink(self, it_xyz, iv_filePREPI, iv_numResidue, iv_moleculeID, iv_moleculeNumber):
        """Check the collisions between the atoms of the new residue and all the existing atoms in the system. If there are no collisions the residue is added to the system.
        
        Input-----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue 
            iv_filePREPI (PREPI): name of the PREPI file for the current residue
            iv_numResidue (int): number of the current residue within the crosslinker molecule
            iv_moleculeID (str): letter representing the name of the crosslinker molecule in the system.
            iv_moleculeNumber (int): index number of the crosslinker molecule considering the whole system as a reference.

        Output----------------------
            collision (bool): Boolean value. True if there were collsions in the system, false if there were not.
        """

# 1-Convert values obtained with function int2car to values in atom PDB format
        it_xyz2      = builder.convertxyz(it_xyz)
        t_aux_atoms  = builder.convertToAtoms(it_xyz2,iv_filePREPI.t_name, iv_filePREPI.name, iv_numResidue, iv_moleculeID, iv_filePREPI.t_topologic, iv_filePREPI.t_numbering,iv_filePREPI.t_con)
        t_aux_atoms2 = copy.deepcopy (t_aux_atoms)
# 2-Check if atoms are inside or outside box and calculate Periodic Boundary Condition if it's necessary
        for natoms in range(0, len(t_aux_atoms)):
            t_aux_atoms[natoms].x, t_aux_atoms[natoms].y, t_aux_atoms[natoms].z = \
                               builder.periodicBoundaryCondition(self.system_pbc.box, t_aux_atoms[natoms].x, \
                                                                                t_aux_atoms[natoms].y, \
                                                                                t_aux_atoms[natoms].z )           
# 3-Check if any atom of new residue makes a collision with existing atoms
        for natoms in t_aux_atoms:    
                if natoms.name != 'DUMM':
                    if(self.coll.checkCollisionD(natoms)==True): 
                        collision_flag = True 
                        return collision_flag        
        collision_flag = False 
# 4-Create residue with atoms and other variables
        v_residue  = residue(t_aux_atoms, iv_filePREPI.name, iv_numResidue, iv_moleculeID)
        v_residue2 = residue(t_aux_atoms2, iv_filePREPI.name, iv_numResidue, iv_moleculeID)
# 5-Delete DUMM lines from atom list (Don't want this lines in file that will be created later)

        v_residue.t_atom  = atom.deleteDUMMlines(v_residue.t_atom)
        for elem in v_residue.t_atom:
            self.coll.addAtomToDictionary(elem)
        v_residue2.t_atom = atom.deleteDUMMlines(v_residue2.t_atom)
        
# 6-Add residue to current system
        self.system_pbc.t_molecules[iv_moleculeNumber].addResidue(v_residue)
        self.system_nopbc.t_molecules[iv_moleculeNumber].addResidue(v_residue2)
        
        return collision_flag
######################################################################################################################
######################################################################################################################
    def __MeasureMolecule(self, it_molecule, iv_moleculeNumber, iv_moleculeID, to_filePREPI):
        """Measures the maximum length of a molecule from the first no-dummy-M atom to the last no-dummy-M atom.
        
        Input-----------------------
            it_molecule (list): list of residues in the molecule. 
            iv_moleculeNumber (int): index number of the molecule in the system.
            iv_moleculeID (str): label of the molecule
            to_filePREPI (list): list of prepi files

        Output----------------------
            distance (float): size of the crosslinker
            growth (list): growth vector of the crosslinker
            it_residues (list): list with all the residues of the crosslinker (each residue has atom objects)
        """
        t_xyz              = [] #list of coordinates of every atom in a residue inculding dummies???
        t_last3m           = [] #list of coordinates of last 3 M-atoms including dummies if necesary
        t_aux_last3m       = [] 
        t_aux_xyz          = []
        t_history_last3m   = [0]
        t_history_xyz      = [0]
        iv_moleculeNumber  = 0 #The first atom of the molecule to measure will be located at the origin.
# ---------------------------------------------------------------------------------------------------------------------
        it_residues = []
        vNumResidue = 0
        
        while vNumResidue < len(it_molecule):
            residuRetries = 1
            t_aux_last3m = copy.deepcopy(t_history_last3m[vNumResidue])
            t_aux_xyz    = copy.deepcopy(t_history_xyz)
            t_xyz, v_filePREPI, t_last3m = builder.__processResidue(self, vNumResidue, it_molecule[vNumResidue], \
                                                                      residuRetries, iv_moleculeNumber, \
                                                                      to_filePREPI, t_aux_xyz, t_aux_last3m)

        #After processResidue, a list with the coordinates of all atoms of the residue is recovered, but no resiude has created yet
            it_xyz2      = builder.convertxyz(t_xyz)
            t_aux_atoms  = builder.convertToAtoms(it_xyz2,v_filePREPI.t_name, v_filePREPI.name, vNumResidue, iv_moleculeID, v_filePREPI.t_topologic, v_filePREPI.t_numbering, v_filePREPI.t_con)
        #Create residue with atoms and other variables
            v_residue  = residue(t_aux_atoms, v_filePREPI.name, vNumResidue, iv_moleculeID)
        #Delete DUMM lines from atom list
            if vNumResidue != 0: #Keep the dummy atoms for the first residue. They will be used as a reference.
                v_residue.t_atom  = atom.deleteDUMMlines(v_residue.t_atom)
            
            it_residues.append(v_residue)
            vNumResidue +=1
            t_history_last3m.append(copy.deepcopy(t_last3m))
            t_history_xyz.append(copy.deepcopy(t_xyz))
        
        res_length = [] #list that will contain the length of each residue of the crosslinker
        res_direction = [] #list that will contain the direction of growth of each crossliker residue taking the first M no dummy atom as a reference

        i=0
        ####Measure the length of every residue       
        for res in it_residues:
            if i == 0:
                pos_M_ini = 3 #The first M atom not dummy of the first residue is the 4th atom
            else:
                pos_M_ini = 0 #Dummy lines were removed for all the other residues-> The first M atom is the first atom of the residue
            
            pos_M_end = self.__getm1atomResidue(res,to_filePREPI) #get the position of the last M-atom of the residue of crosslink  
            if i != 0: pos_M_end = pos_M_end - 3 # from the second to the last residues the 3 first dummy atoms are not included in the list                       
            
            atom_M_ini = res.t_atom[pos_M_ini] #sets the atom object corresponding to the 1st M-atom of the residue of crosslink               
            atom_M_end = res.t_atom[pos_M_end]  #sets the atom object corresponding to the last M-atom of the residue of crosslink
            
            v_atom_M_ini = np.array([atom_M_ini.x, atom_M_ini.y, atom_M_ini.z])
            v_atom_M_end = np.array([atom_M_end.x, atom_M_end.y, atom_M_end.z])
            
            distance = np.linalg.norm(v_atom_M_end-v_atom_M_ini)
            growth = (v_atom_M_end - v_atom_M_ini)/distance
            res_length.append(distance)
            res_direction.append(growth)
            i += 1      
        ####

        ###Measure the length of the crosslink 
        if len(it_residues) > 1:
            res_ini = it_residues[0]
            res_end = it_residues[-1]
            
            pos_M_ini = 3 #The first M atom not dummy of the first residue is the 4th atom           
            pos_M_end = self.__getm1atomResidue(res_end,to_filePREPI) #get the position of the last M-atom of the last residue of crosslink
                        
            atom_M_ini = res_ini.t_atom[pos_M_ini] #sets the atom object corresponding to the 1st M-atom of first residue of crosslink               
            atom_M_end = res_end.t_atom[pos_M_end-3]  #sets the atom object corresponding to the last M-atom of last residue of crosslink, t_atom does not include the 3 first dummy atoms                       
            
            v_atom_M_ini = np.array([atom_M_ini.x, atom_M_ini.y, atom_M_ini.z])
            v_atom_M_end = np.array([atom_M_end.x, atom_M_end.y, atom_M_end.z])
            distance = np.linalg.norm(v_atom_M_end-v_atom_M_ini)
            growth = (v_atom_M_end - v_atom_M_ini)/distance
        else:    
            distance = res_length[0]
            growth = res_direction[0]       

        return distance, growth, it_residues
######################################################################################################################
######################################################################################################################
    def __get_atom_ini_idx(self, iv_molecule_idx, iv_residue_idx, starting_atom_name):
        """Gives the index of starting-type atom where the crosslink will bond and start its growth 

        Input-----------------------
        iv_molecule_idx (int): number representing the index of the random molecule in the system
        iv_residue_idx (int): number representing the index of a random residue in the random molecule
        starting_atom_name (str): string representing the name of the atom where the crosslink will bond and start its growth 

        Output----------------------
        iv_atom_idx (int): number representing the index of the atom where the crosslink will bond and start its growth
        """
        it_atoms = [] #list of atoms where the crosslink will start to grow
            
        atoms_in_residue = self.system_nopbc.t_molecules[iv_molecule_idx].t_residues[iv_residue_idx].t_atom #list with atom objects of the residue

        for i_atom in atoms_in_residue:
            if i_atom.name == starting_atom_name:
                it_atoms.append(atoms_in_residue.index(i_atom))
        if len(it_atoms) == 0: 
            return None
        #RANDOMSEED
        #random.seed(1)
        iv_atom_idx = random.choice(it_atoms)

        return iv_atom_idx
######################################################################################################################
######################################################################################################################
    def __get_atom_end_idx(self, it_mol_res, iv_crosslinkLength, iv_bond_length_CL_start, iv_bond_length_CL_close, atom_ini, mol_ini, closing_atom_name, remeasure_flag):
        """Gives the index of the closing atom where the crosslink finishes its growth located at a certain distance of the starting atom.
        
        Input-----------------------
            it_mol_res (list): list containing the indexes of the [molecules,resiudes] in the system where closing-type atoms are found
            iv_crosslinkLength (float): distance between the first and the last M-atoms of the crosslink molecule in deafault growth
            iv_bond_length_CL_start (float): Bond distance between the first atom of the crosslink and the starting atom
            iv_bond_length_CL_close (float): Bond distance between the last atom of the crosslink and the closing atom
            atom_ini (atom): object of the class atom representing the inital atom to grow the crosslink
            mol_ini (int): index number of the molecule where the starting atom is located
            closing_atom_name (str): string representing the name of the atom where the crosslink will finish its growth i.e. closing atom
            remeasure_flag (bool): Boolean that indicates if the distance will consider all the crosslinks size or only the part of deformation  
        
        Output----------------------
            iv_mol_idx (int): index of the molecule containg the ending hydrogen atom to be substitued to create the crosslink
            iv_res_idx (int): index of the residue containg ending hydrogen atom to be substitued to create the crosslink
            iv_atom_idx (int): index of the ending hydrogen atom to be substitued to create the crosslink
        """
        it_mol_res_atom_end = [] #list of lists with indexes of the molecule, the residue and the atom to be substituted

        for mol_res in it_mol_res:
            mol_idx = mol_res[0]
            res_idx = mol_res[1]

            if mol_idx != mol_ini: #This avoids creating crosslinks between the residues of the same molecule
                #for atom_idx in range(0,len(self.system_nopbc.t_molecules[mol_idx].t_residues[res_idx].t_atom)):
                for atom_idx in range(0,len(self.system_pbc.t_molecules[mol_idx].t_residues[res_idx].t_atom)):    
                    #atom_obj = self.system_nopbc.t_molecules[mol_idx].t_residues[res_idx].t_atom[atom_idx]
                    atom_obj = self.system_pbc.t_molecules[mol_idx].t_residues[res_idx].t_atom[atom_idx]
                    if atom_obj.name == closing_atom_name:
                        #Computes de distance
                        atom_obj_arr = np.array((atom_obj.x,atom_obj.y,atom_obj.z))
                        atom_ini_arr = np.array((atom_ini.x,atom_ini.y,atom_ini.z))
                        distance = np.linalg.norm(atom_obj_arr-atom_ini_arr)
                        if remeasure_flag:
                            #if (distance >= iv_crosslinkLength + iv_bond_length_CL_close) and (distance <= iv_crosslinkLength + iv_bond_length_CL_close*(1+self.system_nopbc.maximum_deviation/100)): #Criteria for searching a new closing atom
                            if (distance >= iv_crosslinkLength + iv_bond_length_CL_close) and (distance <= iv_crosslinkLength + iv_bond_length_CL_close*(1+self.system_pbc.maximum_deviation/100)): #Criteria for searching a new closing atom    
                                it_mol_res_atom_end.append( [mol_idx, res_idx, atom_idx])
                        else:
                            #if (distance >= iv_crosslinkLength + iv_bond_length_CL_start*(1-self.system_nopbc.maximum_deviation/100) + iv_bond_length_CL_close ) and (distance <= iv_crosslinkLength + iv_bond_length_CL_start + (iv_bond_length_CL_close)*(1+self.system_nopbc.maximum_deviation/100)): #Initial criteria for default growth crosslink
                            if (distance >= iv_crosslinkLength + iv_bond_length_CL_start*(1-self.system_pbc.maximum_deviation/100) + iv_bond_length_CL_close ) and (distance <= iv_crosslinkLength + iv_bond_length_CL_start + (iv_bond_length_CL_close)*(1+self.system_pbc.maximum_deviation/100)): #Initial criteria for default growth crosslink
                                it_mol_res_atom_end.append( [mol_idx, res_idx, atom_idx])           
        
        if len(it_mol_res_atom_end) == 0: 
                return None, None, None                           
        #RANDOMSEED
        #random.seed(1)
        rnd_mol_res_atom = random.choice(it_mol_res_atom_end)
        iv_mol_idx = rnd_mol_res_atom[0]
        iv_res_idx = rnd_mol_res_atom[1]
        iv_atom_idx = rnd_mol_res_atom[2]

        return iv_mol_idx, iv_res_idx, iv_atom_idx
######################################################################################################################
######################################################################################################################
    def __Rot_Trans_Residue(self, t_xyz, axis, angle,ref_atom, it_trans_atom, direction = None, initial = False, v_bond_length_CL = 0):
        """Rotates the atoms of residue around an arbitrary axis and translates them to a new reference point.
        
        Input-----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue 
            axis (numpy.ndarray): array representing the orthonormal vector. (between crosslink growth direction and starting to closing (end-ini) atoms direction).
            angle (float): angle in radians between the two vectors.  (between crosslink growth direction and starting to closing (end-ini) atoms direction).
            ref_atom (int): index representing the reference atom in the first residue. 
            it_trans_atom (list): reference coordinates of the atom (starting). The first crosslink atom will be located at a distnce of v_bond_length_CL of this atom in the direection defined by direction
            direction (list): vector direction between the starting and closing atoms
            initial (bool): Whether it is the first residue of the crosslink. 
            v_bond_length_CL (float): Bond distance between the first atom of the crosslink and the starting atom.
        
        Output----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue rotated and translated
        """
        #Recover the original coordinates from the coordinates list and reshape into a matrix of dimension Nx3
        it_xyz = copy.deepcopy(t_xyz)
        
        coor_ini = []
        for element in it_xyz:
            x_ini = element[0]
            y_ini = element[1]
            z_ini = element[2]
            coor_ini.append([x_ini,y_ini,z_ini])     
        coor_ini = np.array(coor_ini).reshape((len(it_xyz),3)) #original coordinates of all residue atoms reshaped into a matrix of dimension Nx3
        #Rotation around the orthonormal axis
        coor_rtn = [self.__rot_arb(p,axis,angle) for p in coor_ini] #atoms rotated coordinates
        
        #List with the rotated coordinates of all the atoms
        for ii in range(len(it_xyz)):
                it_xyz[ii][0] = coor_rtn[ii][0]
                it_xyz[ii][1] = coor_rtn[ii][1]
                it_xyz[ii][2] = coor_rtn[ii][2]
 
        #Diference between the position of every atom of the first residue and the reference/first atom of the residue
        it_delta = []
        for element in it_xyz:
            delta_x = element[0] - it_xyz[ref_atom][0]
            delta_y = element[1] - it_xyz[ref_atom][1]
            delta_z = element[2] - it_xyz[ref_atom][2]
            it_delta.append([delta_x,delta_y,delta_z])              
        
        #######################Translation of the first residue################
        if initial:
            dd = v_bond_length_CL
            x_trans = it_trans_atom[0] + dd*direction[0]
            y_trans = it_trans_atom[1] + dd*direction[1]
            z_trans = it_trans_atom[2] + dd*direction[2]

            it_trans_atom = [x_trans, y_trans, z_trans]

        ################################END####################################

        #Traslation to the new reference point 
        ii = 0
        aux_it_xyz = it_xyz
        it_xyz = []
        for element in aux_it_xyz:
            element[0] = (it_trans_atom[0] + it_delta[ii][0])
            element[1] = (it_trans_atom[1] + it_delta[ii][1])
            element[2] = (it_trans_atom[2] + it_delta[ii][2])
            ii += 1
            pos = np.array([element[0],element[1],element[2]])
            it_xyz.append(pos)
        
        return it_xyz
######################################################################################################################
######################################################################################################################  
    def __ReMeasureMolecule(self, it_molecule, iv_NumResidue, iv_Prepi_res, to_filePREPI, it_last3m, it_xyz,remeasure,counting,dense_counting,collision_counting):
        """Measures the maximum length of a molecule from the first no-dummy-M atom to the last no-dummy-M atom 
        
        Input-----------------------
            it_molecule (list): list of residues in the molecule. 
            iv_NumResidue (int): number of the current residue within the crosslinker molecule.
            iv_Prepi_res (PREPI): name of the PREPI file for the current residue.
            to_filePREPI (list): list of prepi files.
            it_last3m (list) : list of the coordinates of the last 3-M atoms of the current residue.
            it_xyz (list): list with the coordinates of all the atoms of the current residue.
            remeasure (boolean): report the remeasuring parameters.
            counting (boolean): report the number of particles sharing boxes with the chain.
            dense_counting (boolean): report the number of particles sharing boxes with the chain incresing the weight of particles close to the first residues.
            collision_counting (boolean): report the number of potential collisions.
        Output----------------------
            distance (float): size of the crosslinker.
            atom_lastM_ini (atom): atom object corresponding to the last M-atom of first residue of crosslink.
            atom_M_ini (atom) atom object corresponding to the 1st M-atom of first residue of crosslink.
            growth (list): growth vector of the crosslinker.
            
        """
        t_xyz              = [] #list of coordinates of every atom in a residue inculding dummies???
        t_last3m           = [] #list of coordinates of last 3 M-atoms including dummies if necesary
        t_aux_last3m       = [] 
        t_aux_xyz          = []
        #t_history_last3m   = [0]
        t_history_xyz      = [0]
        iv_moleculeNumber  = 999 #The molecule number is not relevant in this routine while it is not 0. 
        iv_moleculeID      = "S" #The ID of the molecule is not relevant in this routine 

        it_residues = []
        vNumResidue = iv_NumResidue
        v_filePREPI = iv_Prepi_res
        t_history_last3m   = []

        reg_count = None
        dense_count = None
        collision_count = None

        for _ in range(vNumResidue + 1):
            t_history_last3m.append([0])

        #Create a residue with the current coordinates of the atoms of reference
        it_xyz2      = builder.convertxyz(it_xyz)
        t_aux_atoms  = builder.convertToAtoms(it_xyz2,v_filePREPI.t_name, v_filePREPI.name, vNumResidue, iv_moleculeID, v_filePREPI.t_topologic, v_filePREPI.t_numbering, v_filePREPI.t_con)
        
        weight_factor = 1000
        if counting:
            t_aux_atoms_c = copy.deepcopy(t_aux_atoms)
            #Place atoms within PBC conditions
            for natoms in range(0, len(t_aux_atoms_c)):
                t_aux_atoms_c[natoms].x, t_aux_atoms_c[natoms].y, t_aux_atoms_c[natoms].z = builder.periodicBoundaryCondition(self.system_pbc.box, t_aux_atoms_c[natoms].x,t_aux_atoms_c[natoms].y,t_aux_atoms_c[natoms].z )   
            
            
            if collision_counting:
                collision_count = 0
                for natoms in t_aux_atoms_c:    
                    if natoms.name != 'DUMM':
                        if(self.coll.checkCollisionD(natoms)==True):
                            collision_count +=1



            reg_count = 0
            dense_count = 0
            for natoms in t_aux_atoms_c:    
                if natoms.name != 'DUMM':
                    reg_count += self.coll.getLengthOfList(natoms)
                    if dense_counting: 
                        dense_count += self.coll.getLengthOfList(natoms)*weight_factor
        
        weight_factor = 100

        
        #Create residue with atoms and other variables
        v_residue  = residue(t_aux_atoms, v_filePREPI.name, vNumResidue, iv_moleculeID)
        #Delete DUMM lines from atom list
        v_residue.t_atom  = atom.deleteDUMMlines(v_residue.t_atom)
        it_residues.append(v_residue)
        #vNumResidue +=1 ####################################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
        t_history_last3m.append(copy.deepcopy(it_last3m))
        #t_history_xyz.append(copy.deepcopy(t_xyz))

        vNumResidue += 1
        while vNumResidue < len(it_molecule):
            residuRetries = 1
            t_aux_last3m = copy.deepcopy(t_history_last3m[vNumResidue])   
            t_aux_xyz    = copy.deepcopy(t_history_xyz)
            t_xyz, v_filePREPI, t_last3m = builder.__processResidue(self, vNumResidue, it_molecule[vNumResidue].name, \
                                                                      residuRetries, iv_moleculeNumber, \
                                                                      to_filePREPI, t_aux_xyz, t_aux_last3m)
        #After processResidue, a list with the coordinates of all atoms of the residue is recovered, but no resiude has created yet
            it_xyz2      = builder.convertxyz(t_xyz)
            t_aux_atoms  = builder.convertToAtoms(it_xyz2,v_filePREPI.t_name, v_filePREPI.name, vNumResidue, iv_moleculeID, v_filePREPI.t_topologic, v_filePREPI.t_numbering, v_filePREPI.t_con)
            
            if counting:
                t_aux_atoms_c = copy.deepcopy(t_aux_atoms)
                #Place atoms within PBC conditions
                for natoms in range(0, len(t_aux_atoms_c)):
                    t_aux_atoms_c[natoms].x, t_aux_atoms_c[natoms].y, t_aux_atoms_c[natoms].z = builder.periodicBoundaryCondition(self.system_pbc.box, t_aux_atoms_c[natoms].x,t_aux_atoms_c[natoms].y,t_aux_atoms_c[natoms].z )   
                
                if collision_counting:
                    for natoms in t_aux_atoms_c:    
                        if natoms.name != 'DUMM':
                            if(self.coll.checkCollisionD(natoms)==True):
                                collision_count +=1
                
                for natoms in t_aux_atoms_c:    
                    if natoms.name != 'DUMM':
                        reg_count += self.coll.getLengthOfList(natoms)
                        if dense_counting: 
                            dense_count += self.coll.getLengthOfList(natoms)*weight_factor
                  
        #Create residue with atoms and other variables
            v_residue  = residue(t_aux_atoms, v_filePREPI.name, vNumResidue, iv_moleculeID)
        #Delete DUMM lines from atom list
            v_residue.t_atom  = atom.deleteDUMMlines(v_residue.t_atom)
            
            it_residues.append(v_residue)
            vNumResidue +=1
            t_history_last3m.append(copy.deepcopy(t_last3m))
            t_history_xyz.append(copy.deepcopy(t_xyz))

            if vNumResidue == 1: 
                weight_factor = 100
            elif vNumResidue == 2: 
                weight_factor = 10
            elif vNumResidue > 2:
                weight_factor = 1

        distance = None
        atom_lastM_ini = None
        atom_M_ini = None
        growth = None


        if remeasure:
            res_ini = it_residues[0]  
            res_end = it_residues[-1]

            pos_M_ini = 0 #Dummy lines were removed for all the residues-> The first M atom is the first atom of the residue
            atom_M_ini = res_ini.t_atom[pos_M_ini] #sets the atom object corresponding to the 1st M-atom of first residue of crosslink               
            
            pos_lastM_ini = self.__getm1atomResidue(res_ini,to_filePREPI) #get the position of the last M-atom of the first residue of crosslink
            atom_lastM_ini = res_ini.t_atom[pos_lastM_ini-3]  #sets the atom object corresponding to the last M-atom of first residue of crosslink, t_atom does not include the 3 first dummy atoms                       

            pos_M_end = self.__getm1atomResidue(res_end,to_filePREPI) #get the position of the last M-atom of the last residue of crosslink
            atom_M_end = res_end.t_atom[pos_M_end-3]  #sets the atom object corresponding to the last M-atom of last residue of crosslink, t_atom does not include the 3 first dummy atoms                       
                
            v_atom_M_ini = np.array([atom_M_ini.x, atom_M_ini.y, atom_M_ini.z])
            v_atom_M_end = np.array([atom_M_end.x, atom_M_end.y, atom_M_end.z])
            distance = np.linalg.norm(v_atom_M_end-v_atom_M_ini)
            growth = (v_atom_M_end - v_atom_M_ini)/distance

        return distance, atom_lastM_ini, atom_M_ini, growth,reg_count,dense_count,collision_count
######################################################################################################################
######################################################################################################################
    def __find_connected_atoms(self, t_atoms , atom_name ):
        """Find all atoms connected to a given atom, directly or inderectly, considering a prepi file as a directed graph.
        Breadth-First Search / Depth-First Search

        Input-----------------------
        t_atoms (list): list of atom objects
        atom_name (str): name of the reference atom to analyse the connectivity
        
        Output----------------------
        t_idx_connecte (list): list of indexes representing atoms connected to the given atom
        start_vertex (int): internal prepi numbering of the given atom
        """
        #Defining a graph (dictionary) with the internal numbering of prepi files as keys (vertices)
        graph={0:[], 1:[],2:[],3:[]} #It is necessary to define the vertices for Dummy atoms
        

        for atom in t_atoms:    
            graph[atom.internal_numbering] = []
            if atom.name == atom_name:
                start_vertex = atom.internal_numbering

        #Add the connectivity information
        for atom in t_atoms:
            t_connectivity = graph[atom.connectivity[0]]
            t_connectivity.append(atom.internal_numbering)
        
        #Breadth-First Search algorithm
        
        visited = set()
        queue = [start_vertex]
        connected_vertices = []
        
        while queue:
            current_vertex = queue.pop(0)
            if current_vertex not in visited:
                visited.add(current_vertex)
                if current_vertex != start_vertex:
                    connected_vertices.append(current_vertex)
                for neighbor in graph.get(current_vertex,[]):
                    if neighbor not in visited:
                        queue.append(neighbor)

        #Express connected vertices values as atom indexes in the atom object list (Not considering the 3 first DUMMY atoms and starting from 0)
        t_idx_connected = []
        for vertix in connected_vertices:
            t_idx_connected.append(int(vertix-4))     


        return t_idx_connected, start_vertex                
######################################################################################################################
######################################################################################################################
    def __getLEaPindex(self,rnd_mol_ini,rnd_res_ini,t_molecules,new_mol_end,new_res_end,crosslinkNumber,t_crosslinks,dendron):
        """Retruns the index in LEaP format of the starting and closing residues and the crosslink.
        Input-----------------------
        rnd_mol_ini (int): number representing the index of the random molecule in the system where the starting atom is located.
        rnd_res_ini (int): number representing the index of the random residue in the system where the starting atom is located.
        t_molecules (list): list of residues.
        new_mol_end (int): number representing the index of the random molecule in the system where the closing atom is located.
        new_res_end (int): number representing the index of the random residue in the system where the closing atom is located.
        crosslinkNumber (int): number representing the index of the crosslink in the system.
        t_crosslink (list): list with residue objects of the crosslinker molecule.
        dendron (boolean): wheter the index is for a dendron chain.
       
        Output----------------------
        res_ini_idx_tleap (int): index of the residue where the starting atom is located in LEaP indexing format.
        res_end_idx_tleap (int): index of the residue where the closing atom is located in LEaP indexing format.
        CL_ini_idx_tleap  (int): index of the first residue of the crosslink in LEaP indexing format.
        CL_end_idx_tleap  (int): index of the last residue of the crosslink in LEaP indexing format.
        """
        
        #Get the index of the starting residue in tleap indexing format

        it_molecules = copy.deepcopy(t_molecules)

        if dendron:
            for elm in t_crosslinks[0:crosslinkNumber]:
                it_molecules.append(elm)
            it_molecules.append(t_crosslinks[crosslinkNumber])

        if rnd_mol_ini != 0:
            res_in_system = 0
            for idx_mol in range (0,rnd_mol_ini):
                res_in_system += len(it_molecules[idx_mol])
            res_ini_idx_tleap = res_in_system + rnd_res_ini + 1 #tleap counts from 1 to N
        else:
            res_ini_idx_tleap = rnd_res_ini + 1 #tleap counts from 1 to N
                                
        #Get the index of the closing residue in tleap indexing format
        if new_mol_end == None and new_res_end == None:
            res_end_idx_tleap = None
        else:    
            if new_mol_end != 0:
                res_in_system = 0
                for idx_mol in range (0,new_mol_end):
                    res_in_system += len(it_molecules[idx_mol])
                res_end_idx_tleap = res_in_system + new_res_end + 1 #tleap counts from 1 to N
            else:
                res_end_idx_tleap = new_res_end + 1 #tleap counts from 1 to N
                                
        #Get the index of the starting Crosslink residue in tleap indexing format
        res_in_mol = 0
        for mol_n in t_molecules:
            res_in_mol += len(mol_n)
                                
        if crosslinkNumber != 0:
            res_in_system = 0
            for idx_mol in range (0,crosslinkNumber):
                res_in_system += len(t_crosslinks[idx_mol])
            CL_ini_idx_tleap = res_in_mol + res_in_system + 1 #tleap counts from 1 to N
            CL_end_idx_tleap = res_in_mol + res_in_system + len(t_crosslinks[crosslinkNumber])  #tleap counts from 1 to N
        else:
            CL_ini_idx_tleap = res_in_mol + 1 #tleap counts from 1 to N                                
            CL_end_idx_tleap = res_in_mol + len(t_crosslinks[crosslinkNumber])  #tleap counts from 1 to N

        return res_ini_idx_tleap, res_end_idx_tleap, CL_ini_idx_tleap, CL_end_idx_tleap        
    
######################################################################################################################
######################################################################################################################
    def __addNewGraft(self, it_molecule, iv_moleculeNumber, iv_moleculeID, to_filePREPI, iv_mol_idx, iv_res_idx, iv_atom_idx, iv_vgrowth_start, iv_growth_cl,atom_ini,t_res_cl,iv_bond_length_CL_start): 
        """
        Adds a new molecule in the system.
        
        Input-----------------------
            it_molecule (list): list of residues in the molecule.
            iv_moleculeNumber (int): index number of the molecule considering the whole system as a reference.
            iv_moleculeID (str): letter representing the name of the molecule in the system.
            to_filePREPI (list): list of all the prepi files.
            iv_mol_idx  (int): number representing the index of the random molecule in the system where the starting atom is located.
            iv_res_idx  (int): number representing the index of the random residue in the system where the starting atom is located.
            iv_atom_idx (int): number representing the index of the starting atom.
            iv_CL_dist (float): Bond distance between the first atom of the crosslink and the starting atom.
            iv_vgrowth_start (numpy.ndarray): array with the unitary growth vector where the starting atom is located.
            iv_growth_cl (numpy.ndarray): array with the unitary growth vector of the crosslink

        Output----------------------
        v_error (int): flag variable showing if the graft chain was created.
        """


        t_xyz              = [] #list of coordinates of every atom in a residue inculding dummies
        v_error            = 0  #Value of error flag variable 
        t_last3m           = [] #list of coordinates of last 3 M-atoms including dummies if necesary
        t_aux_last3m       = [] 
        t_aux_xyz          = []
        t_history_last3m   = [0]
        t_history_xyz      = [0]
        t_residuRetries    = []  #Number of retries of residue number
        collision          = False
# 1-Create new molecule with molecule name from file. Here the molecule object is added with no processed residues
        v_molecule = molecule(iv_moleculeID, iv_moleculeID)
        self.system_pbc.addMolecule(v_molecule)
        v_molecule2 = molecule(iv_moleculeID, iv_moleculeID)
        self.system_nopbc.addMolecule(v_molecule2)
# ---------------------------------------------------------------------------------------------------------------------
        for fillResiduRetries in range(0,len(it_molecule)): t_residuRetries.append(0)
        vNumResidue = 0
        while vNumResidue < len(it_molecule):
            t_residuRetries[vNumResidue] += 1
            if t_residuRetries[vNumResidue] < self.system_pbc.max_retries:   #Check if there are still retries for this residue
                t_aux_last3m = copy.deepcopy(t_history_last3m[vNumResidue]) 
                t_aux_xyz    = copy.deepcopy(t_history_xyz) 

                t_xyz, v_filePREPI, t_last3m = builder.__processResidueGraft(self, vNumResidue, it_molecule[vNumResidue],t_residuRetries[vNumResidue],\
                                                                      to_filePREPI, t_aux_xyz, t_aux_last3m, iv_growth_cl, iv_vgrowth_start,\
                                                                        atom_ini,t_res_cl,iv_bond_length_CL_start)
                
                
                #Molecules objects were created and added to the system before __processResidue. 
                #After processResidue, a list with the coordinates of all atoms of the residue is recovered, but no resiude has created yet

                collision = self.__addNewResidueCrossLink(t_xyz, v_filePREPI, vNumResidue, iv_moleculeID, iv_moleculeNumber)
                
                #In __addNewResidue all the atoms objects were created. If no collisions are detected, the resiudes are created and placed in the system in the same function.
                # Then, if this happens the output of funtcion will be False. 
                # Otherwise, residues are neither created nor placed. And the retries will increment for the current residue

                if collision == False:                              #No collision, process next residue
                    vNumResidue +=1
                    self.system_pbc.state[1] = vNumResidue
                    self.system_nopbc.state[1] = self.system_pbc.state[1]
                    self.bar.add(1)                                  #Update progress bar counter
                    t_history_last3m.append(copy.deepcopy(t_last3m)) # [0, [[xm1,ym1,zm1],...,[xm3,ym3,zm3]],...]
                    t_history_xyz.append(copy.deepcopy(t_xyz))   # [0, [[x1,y1,z1],...,[xn,yn,zn]],...]

            else:
                if vNumResidue == 0:                                #No residue could be placed for this molecule. The molecule could not be built
                    return self.system_pbc, self.system_nopbc, 1
                else:
                    self.printLogs( '  Get previous residue', 2 )
                    t_residuRetries[vNumResidue]    = 0              #Reset retries for this residue
                    vNumResidue -=1                                  #Process again previous residue
                    self.bar.substract(1)                            #Update progress bar counter
                    #All information of the previous residue that was already placed in the system is deleted, so that it will be created again
                    del t_history_last3m[-1]
                    del t_history_xyz[-1]
                    self.coll.removeFromDictionary(self.system_pbc.t_molecules[iv_moleculeNumber].t_residues[-1].t_atom)
                    del self.system_pbc.t_molecules[iv_moleculeNumber].t_residues[-1]
                    del self.system_nopbc.t_molecules[iv_moleculeNumber].t_residues[-1]
# ---------------------------------------------------------------------------------------------------------------------
#v_error will be always 0 unless the program enters to the condition of residue num=0 and retries=maxretries        
        return v_error
######################################################################################################################
######################################################################################################################
    def __processResidueGraft(self, iv_numResidue, iv_residueName, iv_numRetry, it_filePREPI, \
                         it_xyz, it_last3m,vector_growth,vector_start,atom_ini, t_res_cl,v_bond_length_CL_start):
        """Get the coordinates of a residue in a molecule generated from the geometry of the prepi file.
        
        Input-----------------------
            iv_numResidue (int): number of the current residue within the crosslinker molecule
            iv_residueName (str): name of the current residue
            iv_numRetry (int): number of retries for the current residue 
            it_filePREPI (list): list of all the prepi files
            it_xyz (list): list of the coordinates of all the residues placed in the system
            it_last3m (list): list of the coordinates of the last 3-M atoms of  the last residue placed for the crosslinker molecule (void if it is the first residue of the molecule)
            vector_growth (list): vector of growth direction of the crosslinker molecule
            vector_end_ini (list): vector direction between the starting and closing atoms
            atom_ini (atom): starting atom. 
            atom_end_idx (list): list with the indexes of the molecule, the residue and the closing atom.
            atom_end_name (str): closing atom name.
            t_res_cl (list): list with residue objects of the crosslinker molecule. The first one has dummy atoms.
            t_mol_res (list): list containing the indexes of the [molecules,resiudes] in the system where closing-type atoms are found
            v_mol_ini (int): index number of the molecule where the starting atom is located
            v_distortion_angle (float): value of the maximum distortion angle a crosslink residue can rotate to avoid a colission and point to a new closing atom
            v_bond_length_CL_start (float): Bond distance between the first atom of the crosslink and the starting atom
            v_bond_length_CL_close (float): Bond distance between the last atom of the crosslink and the closing atom


        Output----------------------
            it_xyz (list): list with the coordinates of all the atoms of the current residue.
            v_filePREPI (PREPI): name of the PREPI file for the current residue.
            it_last3m (list) : list of the coordinates of the last 3-M atoms of the current residue.
            v_error (int): flag variable showing if the residues were placed. 
            change_flag (bool): wheter a new closing atom has been found.

        """
        for x in range(0, len(it_filePREPI)):                     #Get correct filePREPI value from all filePREPI files
            if it_filePREPI[x].name == iv_residueName: v_filePREPI = it_filePREPI[x]

        if iv_numRetry > 1 : #The second and next attempts to place the next residues. Modifies randomly th,ph from original PREPI info
            t_aux_par = builder.getRandomThetaAndPhi(v_filePREPI.t_par, self.system_pbc.maximum_deviation)
        else: #The first attempt to place the residue, from the second residue, uses r,th,ph from prepi File
            t_aux_par = v_filePREPI.t_par #list with r,angle,dihedral for every atom in PREPI including DUMMY

        if iv_numResidue == 0: #The first residue of the molecule will be placed on a perpendicular vector to the growth vector of the starting atom residue at a distince v_bond_length_CL_start from atom_ini.                  
            density_ref = 10
            xyz_reorient = []
            last3m_reorient = []
            #Stablish a new reference system
            axis_x = vector_start
            axis_z = self.__RotateVector(vector_start,vector_growth,np.pi/2)
            axis_pxz = self.__RotateVector(vector_start,vector_growth,np.pi*1/4)
            axis_nxz = self.__RotateVector(vector_start,vector_growth,np.pi*3/4)
            axis_y = np.cross(vector_start,vector_growth)
            axis_y /= np.linalg.norm(axis_y)

            it_res_cl = copy.deepcopy(t_res_cl)
            it_res0_cl = [] #list that will containt the atoms objects of the first residue of crosslink
                        
            for atom_cl in it_res_cl[0].t_atom:
                x_atom = atom_cl.x
                y_atom = atom_cl.y  
                z_atom = atom_cl.z  
                it_res0_cl.append([x_atom,y_atom,z_atom])
 
            #Align graft vector to axis_x
            angle_xz = np.arccos(np.dot(vector_start,vector_growth))
            it_xyz = self.__Rot_Trans_Residue(it_res0_cl, axis_y, angle_xz,3, [atom_ini.x, atom_ini.y, atom_ini.z], axis_x, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
            it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
            it_align_x = [list(vec) for vec in it_xyz]

            rot_angle = [0.0, np.pi*1/4, np.pi*2/4, np.pi*3/4, np.pi*4/4, np.pi*5/4, np.pi*6/4, np.pi*7/4 ]

            #Rotations on plane x-z    
            for angle_xz in rot_angle:
                v_direction = np.cos(angle_xz)*axis_x + np.sin(angle_xz)*axis_z
                it_xyz = self.__Rot_Trans_Residue(it_align_x, axis_y, angle_xz,3, [atom_ini.x, atom_ini.y, atom_ini.z], v_direction, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
                it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
                _,_,_,_,_,_,coll_c = self.__ReMeasureMolecule( t_res_cl, iv_numResidue, v_filePREPI, it_filePREPI, it_last3m, it_xyz,False,True,True,True)
                if coll_c < density_ref:
                    xyz_reorient = []
                    last3m_reorient = []
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
                    density_ref = coll_c
                elif coll_c == density_ref:
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)

            #Rotations on plane x-y    
            for angle_xy in rot_angle:
                v_direction = np.cos(angle_xy)*axis_x + np.sin(angle_xy)*axis_y
                it_xyz = self.__Rot_Trans_Residue(it_align_x, axis_z, angle_xy,3, [atom_ini.x, atom_ini.y, atom_ini.z], v_direction, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
                it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
                _,_,_,_,_,_,coll_c = self.__ReMeasureMolecule( t_res_cl, iv_numResidue, v_filePREPI, it_filePREPI, it_last3m, it_xyz,False,True,True,True)
                if coll_c < density_ref:
                    xyz_reorient = []
                    last3m_reorient = []
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
                    density_ref = coll_c
                elif coll_c == density_ref:
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
                
            #Align graft vector to axis_z
            angle_zx = -np.arccos(np.dot(vector_growth,axis_z))
            it_xyz = self.__Rot_Trans_Residue(it_res0_cl, axis_y, angle_zx,3, [atom_ini.x, atom_ini.y, atom_ini.z], axis_z, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
            it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
            it_align_z = [list(vec) for vec in it_xyz]
            #Rotations on plane y-z    
            for angle_yz in rot_angle:
                v_direction = np.cos(angle_yz)*axis_y + np.sin(angle_yz)*axis_z
                it_xyz = self.__Rot_Trans_Residue(it_align_z, axis_x, angle_yz,3, [atom_ini.x, atom_ini.y, atom_ini.z], v_direction, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
                it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.            
                _,_,_,_,_,_,coll_c = self.__ReMeasureMolecule( t_res_cl, iv_numResidue, v_filePREPI, it_filePREPI, it_last3m, it_xyz,False,True,True,True)
                if coll_c < density_ref:
                    xyz_reorient = []
                    last3m_reorient = []
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
                    density_ref = coll_c
                elif coll_c == density_ref:
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
            #Align graft vector to axis_pxz-y
            angle_pxz = -np.arccos(np.dot(vector_growth,axis_pxz))
            it_xyz = self.__Rot_Trans_Residue(it_res0_cl, axis_y, angle_pxz,3, [atom_ini.x, atom_ini.y, atom_ini.z], axis_pxz, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
            it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
            it_align_pxz = [list(vec) for vec in it_xyz]
            #Rotations on plane pxz-y    
            for angle_pxzy in rot_angle:
                v_direction = np.cos(angle_pxzy)*axis_pxz + np.sin(angle_pxzy)*axis_y
                it_xyz = self.__Rot_Trans_Residue(it_align_pxz, axis_nxz, angle_pxzy,3, [atom_ini.x, atom_ini.y, atom_ini.z], v_direction, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
                it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
                _,_,_,_,_,_,coll_c = self.__ReMeasureMolecule( t_res_cl, iv_numResidue, v_filePREPI, it_filePREPI, it_last3m, it_xyz,False,True,True,True)
                if coll_c < density_ref:
                    xyz_reorient = []
                    last3m_reorient = []
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
                    density_ref = coll_c
                elif coll_c == density_ref:
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
            #Align graft vector to axis_nxz-y
            angle_nxz = -np.arccos(np.dot(vector_growth,axis_nxz))
            it_xyz = self.__Rot_Trans_Residue(it_res0_cl, axis_y, angle_nxz,3, [atom_ini.x, atom_ini.y, atom_ini.z], axis_nxz, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
            it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
            it_align_nxz = [list(vec) for vec in it_xyz]

            #Rotations on plane nxz-y    
            for angle_nxzy in rot_angle:
                v_direction = np.cos(angle_nxzy)*axis_nxz + np.sin(angle_nxzy)*axis_y
                it_xyz = self.__Rot_Trans_Residue(it_align_nxz, axis_pxz, angle_nxzy,3, [atom_ini.x, atom_ini.y, atom_ini.z], v_direction, True, v_bond_length_CL_start )#3 is the first no dummy M atom            
                it_last3m = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, it_last3m) #get the coordinates of last 3M atoms, including DUMMY if necesary.This list is initialized here always.
                _,_,_,_,_,_,coll_c = self.__ReMeasureMolecule( t_res_cl, iv_numResidue, v_filePREPI, it_filePREPI, it_last3m, it_xyz,False,True,True,True)
                if coll_c < density_ref:
                    xyz_reorient = []
                    last3m_reorient = []
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
                    density_ref = coll_c
                elif coll_c == density_ref:
                    xyz_reorient .append(it_xyz)
                    last3m_reorient.append(it_last3m)
            if len(xyz_reorient) == 1:
                it_xyz = xyz_reorient [0]
                it_last3m = last3m_reorient [0] 
            else:
                idx_rnd = random.randrange(len(xyz_reorient))
                it_xyz = xyz_reorient [idx_rnd]
                it_last3m = last3m_reorient [idx_rnd]

        else:  #iv_numResidue > 0
            t_old_last3m = copy.deepcopy (it_last3m) #it_last3m always is initialized for the first residue of every molecule. The second residue will use info of M atoms of first residue if 2nd residue has less than 3 M-atoms
            it_xyz       = int2car.int2car(v_filePREPI.t_con, t_aux_par, it_last3m)#Gets the coordinates using the modified th,ph if necessary. Uses last 3 M-atoms as reference. 
                                                                                   #It places only no dummy atoms, the first 3 in the list are the last 3 M-atoms of the previous residue
                                                                                   #These first 3 M-atoms will be labeled as dummies by the convertToAtoms method
            it_last3m    = builder.getLast3M(v_filePREPI.t_mLines, it_xyz, iv_numResidue, t_old_last3m)#updates the last 3 M-atoms coordintaes, considering the info of the former residue.
                        
        return it_xyz, v_filePREPI, it_last3m
######################################################################################################################
######################################################################################################################
    def __RotateVector(self,v_g,v_s,angle):
        """Rotates an unit vector by an angle on the plane defined by non collinears vectors v_a and v_b.
        
        Input-----------------------
            v_g (numpy.ndarray): unit vector to rotate on the plane defined by v_a and v_b.
            v_s (numpy.ndarray): unit vector that defines a plane with v_a.
            angle (float): rotation angle in radians.

        Output-----------------------
            v_gr (numpy.ndarray): rotated unit vector.
        """

        #Perpendicular component on the plane AB

        if (v_g[0] == v_s[0] and v_g[1] == v_s[1] and v_g[2] == v_s[2]) or (v_g[0] == -v_s[0] and v_g[1] == -v_s[1] and v_g[2] == -v_s[2]):
            v_p = np.array([1.0,0.0,0.0]) if v_g != np.array([1.0,0.0,0.0]) else np.array([0.0,1.0,0.0])
        else:
            v_p = v_s -np.dot(v_s,v_g)*v_g
            v_p /= np.linalg.norm(v_p)

        #Rotate the vector

        v_gr = np.cos(angle)*v_g +np.sin(angle)*v_p

        return v_gr
######################################################################################################################
######################################################################################################################


