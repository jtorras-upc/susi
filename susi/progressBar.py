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

from threading import Timer
import time
import sys

class progressBar:
######################################################################################################################
######################################################################################################################

    def __init__ (self, _max, _prefix, _sufix, _size):
        """Class of object to represent a progress bar while constructing the system

        Input-----------------------
            _max (int): Number of maximum of items to be count.
            _prefix (str): Prefix to be shown before progress bar.
            _sufix (str): Sufix to be shown after progress bar.
            _size (int): Progress bar lenght.
        """
        self.count  = 0              # Items counter
        self.delta  = 0.0            # Delta time in seconds counting
        self.inctime= 0.5            # Tick time in seconds
        self.end    = False          # Counting ended ?
        self.max    = _max           # Maximum number of items to count
        self.prefix = _prefix        # Text before Progress Bar
        self.sufix  = _sufix         # Item units to be count
        self.size   = _size          # Progress bar size
        it          = "."*_size
        #self.t = Timer(self.inctime, self.show)
        self.show()
        sys.stdout.flush()
    
######################################################################################################################
######################################################################################################################
    def show(self):
        """Show current status of the progress bar
        """
        x = int(self.size*self.count/self.max)
        sys.stdout.write("%s[%s%s] %6i/%6i %s [delta %7.1f s]\r" % 
                 (self.prefix, "#"*x, "."*(self.size-x), self.count, self.max, self.sufix,self.delta ))
        sys.stdout.flush()
        self.delta += self.inctime
        if not self.end :      # Start clock tick timer only if counting is not terminated 
            self.t = Timer(self.inctime, self.show)
            self.t.start()
    
######################################################################################################################
######################################################################################################################
    def add(self, num):
        """Increases number of units the bar counter and update view.

        Input-----------------------
        num (int): Increasing value.
        """
        self.count += num

######################################################################################################################
######################################################################################################################
    def substract(self, num):
        """Decreases number of units the bar counter and update view.

        Input-----------------------
        num (int): Decreasing value.
        """
        self.count -= num

######################################################################################################################
######################################################################################################################
    def terminate(self):
        """Terminates the progress bar counter.
        """
        self.end = True
        time.sleep(self.inctime)
        sys.stdout.write("\n### Process completed ###\n")
