#!/usr/bin/env python
###############################################################################
#                                                                             #
#    ParseM.py                                                                #
#                                                                             #
#    Generic parsers for bioinformatics formats                               #
#                                                                             #
#    Copyright (C) Michael Imelfort, Donovan Parks                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort, Donovan Parks"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################
from BamParser import BamParser as BP

class TemplateClass():
    """Utilities wrapper"""
    def __init__(self): pass

    def sayHi(self):
        print("ParseM doesn't so anything on it's own. Perhaps we could put some demos here?" )
        fred = BP()

###############################################################################
###############################################################################
###############################################################################
###############################################################################
