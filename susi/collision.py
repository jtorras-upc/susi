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

class collision:
    __version__ = '2.0.0'
######################################################################################################################
######################################################################################################################
    def __init__(self,box,minimum_distance,npart):
        """Class that transforms atoms coordinates in indexes within of subdivisions of the main box.

        Input-----------------------
        box (list): Box of the system.
        minimum_distance (float): Minimum distance of collision.
        npart (int): The number of partitions along every axis
        """
        self.box                     = box                                #Periodical box where our system is contained
        self.minimum_distance        = minimum_distance                   #Minimum distance of collision
        self.minimum_square_distance = minimum_distance*minimum_distance  #The square of the previous distance
        self.npart                   = npart                              #The number of partitions along every axis
        self.xpart                   = 0                                  #The length of every X partition
        self.ypart                   = 0                                  #The length of every Y partition
        self.zpart                   = 0                                  #The length of every Z partition
        self.LHX                     = 0                                  #Half distance of the box's X dimension
        self.LHY                     = 0                                  #Half distance of the box's Y dimension 	
        self.LHZ                     = 0                                  #Half distance of the box's Z dimension
        self.Dictionary              = {}                                 #The hashtable where we have the list of every 
                                                                          #  atom depending on which part of the box it is 
######################################################################################################################
######################################################################################################################
    def createDic(self):
        """ Generate all keys needed in the dictionary and their empty list to store the atoms in the corresponding partition depending on its posisition. 
        It obtains and assigns the values of the lengths xpart, ypart, zpart, LHX, LHY and LHZ initialized at 0 in the init step.

        Output----------------------
        Dictionary (dict): dictionary to store the atoms in the corresponding partition depending on its posisition.
        """
        self.xpart=self.box[0]/self.npart       #Create and save the correct value for the xpart                      
        self.ypart=self.box[2]/self.npart       #Create and save the correct value for the ypart
        self.zpart=self.box[4]/self.npart       #Create and save the correct value for the zpart
        self.LHX=self.box[0]/2                  #Create and save the correct value for the LHX
        self.LHY=self.box[2]/2                  #Create and save the correct value for the LHY
        self.LHZ=self.box[4]/2                  #Create and save the correct value for the LHZ
        for i in range (self.npart):            #This bucle generate all the keys and empty lists needed, it's 
                                                        #  important to say that the multipliers to i,j,k are not random, they are
                                                        #  choosen with this values to never get the same key with different partitions.
            for j in range (self.npart):
                for k in range(self.npart):
                    self.Dictionary[((self.npart*3)**2)*i+(self.npart*3)*j+k]=[]
        return self.Dictionary
######################################################################################################################
######################################################################################################################
    def checkCollisionD(self,atom):
        """Check the collision of the input atom with the atoms living in the same fragment of space.

        Input-----------------------
        atom (atom): the atom object which one we are cheking it's not collisioning.
        Output----------------------
        flag (bool): wheter a collision is detected.
        """
        b, c, d=collision.getKeyFromAtomPosition(self,atom)                #Obtains the values to generate the corresponding key due to its position
        a= self.Dictionary[((self.npart*3)**2)*b+(self.npart*3)*c+d]       #List of the atoms living in the key box acording to the atom position
        if (collision.checkDistance(self,a,atom) == True):                 #Check if there is collision with the atoms on the previous list and our atom
            return True	                                           #Collision detected = End of the function
        if (collision.nextBox(self,b,c,d,atom) == True):                   #Check if our atom is close enough to the wall box. If it is close enough,
                                                                                   #  it will check collisions on the box which our atom is close to.
            return True                                                #Collision detected = End of the function
        return False                                                       #No collision detected in all checks. The atom is not collisioning = End of the function
######################################################################################################################
######################################################################################################################
    @staticmethod   
    def checkDistance(self,list_of_atoms,atom):
        """Once we choose a key from dictionary cheks if the distance of our atom respect every atom in the list is smaller than the collision distance.
        
        Input-----------------------
        list_of_atoms (list): The list of atom objects of the corresponding key of the dictionary.
        atom (atom): The atom object which we are cheking if collides.
        Output----------------------
        flag (bool): wheter a collision is detected.
        """

        if (len(list_of_atoms)>0):                                     #If the list is empty there is nothing to check
            for elem in list_of_atoms:                             #Loop that goes along every atom of the list
                x=elem.x                                       #X position from the atom of the list 
                y=elem.y                                       #Y position from the atom of the list
                z=elem.z                                       #Z position from the atom of the list
                dist=(x-atom.x)**2+(y-atom.y)**2+(z-atom.z)**2 #Calculation of the square of the distance between two atoms 
                if dist< self.minimum_square_distance:         #Collision condition
                    return True                            #Collision detected = End of the function
        return False                                                   #No collision detected = End of the function
######################################################################################################################
######################################################################################################################
    @staticmethod
    def nextBox(self,b,c,d,atom):
        """Check if exists any collision with atoms of the boxes colliding with current atom box in the case it is close enough to the wall of the box.
        
        Input-----------------------
        b (int): The X identifier of the key.
        c (int): The Y identifier of the key.
        d (int): The Z identifier of the key.
        atom (atom): The atom which we are cheking if collides.
        Output----------------------
        flag (bool): wheter a collision is detected.
        """
        vector1=[b,c,d]                                                                                         #Take the keys and put them in a vector
        vector2=[(atom.x+self.LHX),(atom.y+self.LHY),(atom.z+self.LHZ)]                                         #Take the components of the atom position and displace them to get all time positions between 0 and L, instead of -L/2 and L/2, and put them in a vector
        vector3=[self.xpart,self.ypart,self.zpart]                                                              #Take the partitions of every axis in a vector
        for i in range (len(vector1)):                                                                          #Loop along every dimension
            vector4=vector1[:]                                                                              #Copy the key vector in a new memory place
            if (abs(vector1[i]*vector3[i]-vector2[i])<self.minimum_distance):                               #We check if the distance between the atom and one of the 2 walls of the box for each dimension is near
                if vector1[i]==(0):                                                                     #Check if we are near of the wall with the key 0
                    vector4[i]=self.npart-1                                                         #If it is true the new key is the highest
                else:
                    vector4[i]=vector1[i]-1                                                         #If not, we just decrease 1 unit the key
                a=self.Dictionary[((self.npart*3)**2)*vector4[0]+(self.npart*3)*vector4[1]+vector4[2]]  #We take the list of the atoms from the box which wall we are near to
                if collision.checkDistance(self,a,atom)==True:                                          #We check if there is collision
                    return True                                                                     #Collision detected = End of the function
            if (abs((vector1[i]+1)*vector3[i]-vector2[i])<self.minimum_distance):                           #We check if we are near to the other wall of the box in the same dimension
                if vector1[i]==(self.npart-1):                                                          #If the key is the highest the new box we will check is 0
                    vector4[i]=0
                else:                                                                                   #If not we just add 1 to the key number
                    vector4[i]=vector1[i]+1
                a=self.Dictionary[((self.npart*3)**2)*vector4[0]+(self.npart*3)*vector4[1]+vector4[2]]  #The new list corresponding to the key
                if collision.checkDistance(self,a,atom)==True:                                          #Collision check
                    return True                                                                     #Collision detected = End of the funcition
        return False                                                                                            #No collision detected in all cases = End of the function
#####################################################################################################################
######################################################################################################################
    def addAtomToDictionary(self,atom):
        """Add the atom in the input to its corresponding key in the dictionary.

        Input-----------------------
        atom (atom):  Atom object to be added to the dictionary.
        """
        b, c, d=collision.getKeyFromAtomPosition(self,atom)               #We obtain the key of the atom corresponding to its position                                  
        a= self.Dictionary[((self.npart*3)**2)*b+(self.npart*3)*c+d]      #With the key, we obtain the list of atoms
        a.append(atom)                                                    #We append our atom to this list
######################################################################################################################
######################################################################################################################
    def removeFromDictionary(self,t_atom):
        """Remove all atoms contained in the list(t_atom) from the dictionary.
        
        Input-----------------------
        t_atom (list): List with atom objects.
        """
        for elem in t_atom:                                                       #Bucle for every atom in the list
            b, c, d=collision.getKeyFromAtomPosition(self,elem)               #We obtain the key of the atom corresponding to its position
            a= self.Dictionary[((self.npart*3)**2)*b+(self.npart*3)*c+d]      #With the key, we obtain the list of atoms
            a.remove(elem)                                                    #We remove the atom from the list
######################################################################################################################
######################################################################################################################
    @staticmethod
    def getKeyFromAtomPosition(self,atom):
        """From the atom x,y,z position we obtain the key of the dictionary 
       
        Input-----------------------
        atom (atom):  Atom object we want to get the key.
        Output----------------------
        b (int): Key corresponding to X dimension.
        c (int): Key corresponding to Y dimension.
        d (int): Key corresponding to Z dimension
        """
        b=int((atom.x+self.LHX)/self.xpart)       #Key corresponding to x dimension, we calculate using displace and then dividing with the xpart
        if (b == self.npart):                     #If we are on the highest limit we put the atom in the highest box
            b=self.npart-1                    
        c=int((atom.y+self.LHY)/self.ypart)       #Key corresponding to y dimension, same form to calculate and same check
        if (c == self.npart):                     
            c=self.npart-1                    
        d=int((atom.z+self.LHZ)/self.zpart)       #Key corresponding to z dimension, same form to calculate and same check
        if (d == self.npart):                     
            d=self.npart-1                   
        return b, c, d                            #Return the key composed by the other 3 keys

#######################################################################################################################
#######################################################################################################################
    def getLengthOfList(self,atom):
        """Get the number of atoms in the same sub-box of a specific atom.

        Input-----------------------
        atom (atom):  Atom object we want to obtain the number of atoms in its box.
        Output----------------------
        length (int): Number of atoms in the same sub-box of a specific atom.
        """
        b, c, d=collision.getKeyFromAtomPosition(self,atom)                #Obtains the values to generate the corresponding key due to its position
        a= self.Dictionary[((self.npart*3)**2)*b+(self.npart*3)*c+d]       #With the key, we obtain the list of atoms
        length=len(a)
        return length                            #Return the length of list
#######################################################################################################################
#######################################################################################################################
    def getAtomsOfList(self,atom):
        """Get the names of atoms in the same sub-box of a specific atom.

        Input-----------------------
        atom (atom):  Atom object we want to obtain the number of atoms in its box.
        Output----------------------
        names (list): Names of atoms in the same sub-box of a specific atom.
        """
        b, c, d=collision.getKeyFromAtomPosition(self,atom)                #Obtains the values to generate the corresponding key due to its position
        a= self.Dictionary[((self.npart*3)**2)*b+(self.npart*3)*c+d]       #With the key, we obtain the list of atoms
        
        names = [ element.name for element in a]
        
        return names         
#######################################################################################################################
#######################################################################################################################
    def checkAtomsinList(self,atom_check,atom):
        """Checks if atoms_check is already in the list of atom.

        Input-----------------------
        atom_check (atom): Atom to check whether is in the list
        atom (atom):  Atom object we want to obtain the number of atoms in its box.
        Output----------------------
        names (list): Names of atoms in the same sub-box of a specific atom.
        """
        b, c, d=collision.getKeyFromAtomPosition(self,atom)                #Obtains the values to generate the corresponding key due to its position
        a = self.Dictionary[((self.npart*3)**2)*b+(self.npart*3)*c+d]       #With the key, we obtain the list of atoms
        
        for element in a:
            if element == atom_check: return True
        
        return False                                   
