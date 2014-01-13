#!/usr/bin/env python
###############################################################################
#                                                                             #
#    ContigParser.py                                                          #
#                                                                             #
#    Class for parsing contigs                                                #
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

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readFasta(self, fp): # this is a generator function
        header = None
        seq = None
        while True:
            for l in fp:
                if l[0] == '>': # fasta header line
                    if header is not None:
                        # we have reached a new sequence
                        yield header, "".join(seq)
                    header = l.rstrip()[1:].partition(" ")[0] # save the header we just saw
                    seq = []
                else:
                    seq.append(l.rstrip())
            # anything left in the barrel?
            if header is not None:
                yield header, "".join(seq)
            break

###############################################################################
###############################################################################
###############################################################################
###############################################################################
