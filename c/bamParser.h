//#############################################################################
//
//   bamParser.h
//   
//   Determine average coverage values and linking read pairs
//
//   Copyright (C) Michael Imelfort
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//#############################################################################

#ifndef PM_BAM_PARSER_H
  #define PM_BAM_PARSER_H

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// htslib
#include "htslib/bgzf.h"
#include "htslib/sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"

typedef BGZF bamFile;

#ifdef __cplusplus
extern "C" {
#endif

extern int vomit(int fred);

/*! @typedef
 @abstract Auxiliary data structure used in read_bam
 @field fp the file handler
 @field iter NULL if a region not specified
 @field min_mapQ mapQ filter
 @field min_len length filter
 */
typedef struct {                    // 
    bamFile *fp;                    // the file handler
    hts_itr_t *iter;                // NULL if a region not specified
    int min_mapQ, min_len;          // mapQ filter; length filter
} aux_t;

/*! @typedef
 @abstract Structure for returning mapping results
 @field plp_bp number of bases piled up on each contig
 @field contig_names names of the reference sequences
 @field contig_lengths lengths of the referene sequences
 @field contig_length_correctors corrections to contig lengths used when doing outlier coverage
 @field num_bams number of BAM files parsed
 @field is_links_included are links being calculated
 @field is_outlier_coverage is outlier adjusted coverage being calculated
 @field is_ignore_supps are supplementary alignments being ignored
 @field num_contigs number of reference sequences
 @field links linking pairs
 */
typedef struct {
    uint32_t ** plp_bp;
    uint32_t * contig_lengths;
    uint32_t ** contig_length_correctors;
    uint32_t is_links_included:1, is_outlier_coverage:1, is_ignore_supps:1, num_bams:29;
    uint32_t num_contigs;
    char ** contig_names;
    cfuhash_table_t * links;
} PM_mapping_results;

int read_bam(void *data,
             bam1_t *b);

/*!
 * @abstract Initialise the mapping results struct
 * 
 * @param MR  mapping results struct to initialise
 * @param BAM_header  htslib BAM header
 * @param numBams  number of BAM files to parse
 * @param doOutlierCoverage  set to 1 if should initialise contig_length_correctors
 * @param doLinks  1 if links should be calculated
 * @param ignoreSuppAlignments  only use primary alignments
 * @return void
 * 
 * @discussion If you call this function then you MUST call destroy_MR
 * when you're done.
 */
void init_MR(PM_mapping_results * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             int doLinks,
             int doOutlierCoverage,
             int ignoreSuppAlignments
);

/*!
 * @abstract Free all the memory calloced in init_MR
 * 
 * @param  MR  mapping results struct to destroy
 * @return void
 * 
 * @discussion If you call this function then you MUST call destroy_MR
 * when you're done.
 */
void destroy_MR(PM_mapping_results * MR);

/*!
 * @abstract Initialise the mapping results struct
 * 
 * @param numBams  number of BAM files to parse
 * @param baseQ  base quality threshold
 * @param mapQ  mapping quality threshold
 * @param minLen  min query length
 * @param doLinks  1 if links should be calculated
 * @param ignoreSuppAlignments  only use primary alignments
 * @param doOutlierCoverage  set to 1 if should initialise contig_length_correctors
 * @param bamFiles  filenames of BAM files ro parse
 * @param MR  mapping results struct to write to
 * @return 0 for success
 * 
 * @discussion This function expects MR to be a null pointer. It calls 
 * init_MR and stores info accordingly. TL;DR If you call this function
 * then you MUST call destroy_MR when you're done.
 * 
 */
int parseCoverageAndLinks(int numBams,
                          int baseQ,
                          int mapQ,
                          int minLen, 
                          int doLinks,
                          int ignoreSuppAlignments,
                          int doOutlierCoverage,
                          char* bamFiles[],
                          PM_mapping_results * MR);
                              
/*!
 * @abstract Adjust (reduce) the number of piled-up bases along a contig
 * 
 * @param  MR  mapping results struct to write to
 * @param  position_holder  array of pileup depths
 * @param  tid  contig currently being processed
 * @param  doOutlierCoverage  remove effects fo very high or very low regions
 * @return void
 * 
 * @discussion This function expects MR to be initialised.
 * it can change the values of contig_length_correctors and plp_bp
 */
void adjustPlpBp(PM_mapping_results * MR,
                 uint32_t ** position_holder,
                 int tid);

/*!
 * @abstract Calculate the coverage for each contig for each BAM
 * 
 * @param  MR  mapping results struct with mapping info
 * @return matrix of floats (rows = contigs, cols = BAMs)
 * 
 * @discussion This function expects MR to be initialised.
 * NOTE: YOU are responsible for freeing the return value
 * recommended method is to use destroyCoverages
 */
float ** calculateCoverages(PM_mapping_results * MR);

/*!
 * @abstract Calculate the coverage for each contig for each BAM
 * 
 * @param covs array to destroy
 * @param numContigs number of rows in array
 * @return void
 */
void destroyCoverages(float ** covs, int numContigs);
    
    /***********************
    *** PRINTING AND I/O ***
    ***********************/

void print_MR(PM_mapping_results * MR);

#ifdef __cplusplus
}
#endif

#endif // PM_BAM_PARSER_H