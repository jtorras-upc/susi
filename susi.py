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

import argparse
from susi import * 

#########################################################################
#usage: Mainbuilder.py [-h] [-t] [-v] -f FILE
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -t, --testmode        Activate test mode, file.pdb won't be created
#  -v, --verbosity       Increase output verbosity
#  -f FILE, --file FILE  Input file with system properties
#########################################################################

print(messages.get_text_copyright())

parser = argparse.ArgumentParser()

parser.add_argument("-t","--testmode", action='store_true',
                    help="Activate test mode, file.pdb won't be created")

parser.add_argument("-v", "--verbosity", action="count",
                    help="Increase output verbosity",default=0)

parser.add_argument("-f","--file",
                    help="Input file with system properties",
		    required=True)
parser.add_argument("-r","--restart",
                    help="Restart file",default="")
args = parser.parse_args()

if args.testmode:
    test = "X"
else:
    test = ""

if args.verbosity > 0:
    print("Verbosity level {} activated".format(args.verbosity))
    
#builder(args.file, args.verbosity, test)
builder(args.file,args.verbosity,test,args.restart)
