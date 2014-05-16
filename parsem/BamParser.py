#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamParser.py                                                             #
#                                                                             #
#    Class for parsing BAM files                                              #
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
import os
import ctypes as c
import pkg_resources

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# fields defined in cfuhash.c but not accessed at this level
class cfuhash_table_t(c.Structure):
    pass

# mapping results structure
"""
typedef struct {
    uint32_t ** plp_bp;
    uint32_t * contig_lengths;
    uint32_t ** contig_length_correctors;
    uint32_t num_bams;
    uint32_t num_contigs;
    char ** contig_names;
    char ** bam_file_names;
    int is_links_included;
    int is_outlier_coverage;
    int is_ignore_supps;
    cfuhash_table_t * links;
} PM_mapping_results;
"""
class PM_mapping_results(c.Structure):
    _fields_ = [("plp_bp", c.POINTER(c.POINTER(c.c_uint32))),
                ("contig_lengths",c.POINTER(c.c_uint32)),
                ("contig_length_correctors",c.POINTER(c.POINTER(c.c_uint32))),
                ("num_bams",c.c_uint32),
                ("num_contigs",c.c_uint32),
                ("contig_names",c.POINTER(c.POINTER(c.c_char))),
                ("bam_file_names",c.POINTER(c.POINTER(c.c_char))),
                ("is_links_included",c.c_int),
                ("is_outlier_coverage",c.c_int),
                ("is_ignore_supps",c.c_int),
                ("links",c.POINTER(cfuhash_table_t))
                ]

class BamParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self):
        # load the c library
        package_dir, filename = os.path.split(__file__)
        package_dir = os.path.abspath(package_dir)
        package_dir = package_dir.replace("parsem","" )
        c_lib = os.path.join(package_dir, 'c', 'bam', 'libPMBam.a')
        self.libPMBam = c.cdll.LoadLibrary(c_lib)

        #---------------------------------
        # import C functions
        #---------------------------------
        self.merge_MR = self.libPMBam.merge_MRs
        """
        @abstract Merge the contents of MR_B into MR_A

        @param  MR_A  mapping results struct to copy to
        @param  MR_B  mapping results struct to copy from
        @return void
        @discussion MR_B remains unchanged.
        MR_A is updated to include all the info contained in MR_B

        void merge_MRs(PM_mapping_results * MR_A, PM_mapping_results * MR_B);
        """

        self.destroy_MR = self.libPMBam.destroy_MR
        """
        @abstract Free all the memory calloced in init_MR

        @param  MR  mapping results struct to destroy
        @return void

        void destroy_MR(PM_mapping_results * MR)
        """

        self.parseCoverageAndLinks = self.libPMBam.parseCoverageAndLinks
        """
        @abstract Initialise the mapping results struct <- read in the BAM files

        @param numBams  number of BAM files to parse
        @param baseQ  base quality threshold
        @param mapQ  mapping quality threshold
        @param minLen  min query length
        @param doLinks  1 if links should be calculated
        @param ignoreSuppAlignments  only use primary alignments
        @param doOutlierCoverage  set to 1 if should initialise contig_length_correctors
        @param bamFiles  filenames of BAM files to parse
        @param MR  mapping results struct to write to
        @return 0 for success

        @discussion This function expects MR to be a null pointer. It calls
        init_MR and stores info accordingly. TL;DR If you call this function
        then you MUST call destroy_MR when you're done.

        int parseCoverageAndLinks(int numBams,
                                  int baseQ,
                                  int mapQ,
                                  int minLen,
                                  int doLinks,
                                  int ignoreSuppAlignments,
                                  int doOutlierCoverage,
                                  char* bamFiles[],
                                  PM_mapping_results * MR
                                 )
        """

        self.adjustPlpBp = self.libPMBam.adjustPlpBp
        """
        @abstract Adjust (reduce) the number of piled-up bases along a contig

        @param  MR  mapping results struct to write to
        @param  position_holder  array of pileup depths
        @param  tid  contig currently being processed
        @param  doOutlierCoverage  remove effects fo very high or very low regions
        @return void

        @discussion This function expects MR to be initialised.
        it can change the values of contig_length_correctors and plp_bp

        void adjustPlpBp(PM_mapping_results * MR,
                         uint32_t ** position_holder,
                         int tid)
        """

        self.calculateCoverages = self.libPMBam.calculateCoverages
        """
        @abstract Calculate the coverage for each contig for each BAM

        @param  MR  mapping results struct with mapping info
        @return matrix of floats (rows = contigs, cols = BAMs)

        @discussion This function expects MR to be initialised.
        NOTE: YOU are responsible for freeing the return value
        recommended method is to use destroyCoverages

        float ** calculateCoverages(PM_mapping_results * MR);
        """

        self.destroyCoverages = self.libPMBam.destroyCoverages
        """
        @abstract Destroy the coverages structure

        @param covs array to destroy
        @param numContigs number of rows in array
        @return void

        void destroyCoverages(float ** covs, int numContigs)
        """

        self.print_MR = self.libPMBam.print_MR
        """
        @abstract Print the contents of the MR struct

        @param  MR   mapping results struct with mapping info

        void print_MR(PM_mapping_results * MR)
        """

###############################################################################
###############################################################################
###############################################################################
###############################################################################
