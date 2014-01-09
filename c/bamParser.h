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

// htslib
#include "htslib/bgzf.h"
#include "htslib/sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"

typedef BGZF bamFile;

/*! @typedef
 @ *a*bstract Auxiliary data structure used in read_bam
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
 @ a*bstract Structure for returning mapping results
 @field total_reads number of reads mapped to each contig
 @field contig_names names of the reference sequences
 @field contig_lengths lengths of the referene sequences
 @field num_bams number of BAM files parsed
 @field num_contigs number of reference sequences
 @field links linking pairs
 */
typedef struct {
    int ** total_reads;
    char ** contig_names;
    float * contig_lengths;
    int num_bams;
    int num_contigs;
    cfuhash_table_t * links;
} PM_mapping_results;

#ifdef __cplusplus
extern "C" {
#endif

int read_bam(void *data,
             bam1_t *b);

/*!
 * @abstract Initialise the mapping results struct
 * 
 * @param  MR  mapping results struct to initialise
 * @param  BAM_header  htslib BAM header
 * @param  numBams  number of BAM files we'll be parsing
 * @return void
 * 
 * @discussion If you call this function then you MUST call destroy_MR
 * when you're done.
 */
void init_MR(PM_mapping_results * MR,
             bam_hdr_t * BAM_header,
             int numBams);

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
 * @param  numBams  number of BAM files we'll be parsing
 * @param  baseQ  base quality threshold
 * @param  mapQ  mapping quality threshold
 * @param  minLen  min query length
 * @param  doLinks  number of BAM files we'll be parsing
 * @param  ignoreSuppAlignments  only use primary alignments
 * @param  doOutlierCoverage  remove effects fo very high or very low regions
 * @param  bamFiles  filenames of BAM files we'll be parsing
 * @param  MR  mapping results struct to write to
 * @return 0 for success
 * 
 */
int processBams(int numBams,
                int baseQ,
                int mapQ,
                int minLen, 
                int doLinks,
                int ignoreSuppAlignments,
                int doOutlierCoverage,
                char* bamFiles[],
                PM_mapping_results * MR);

    /***********************
    *** PRINTING AND I/O ***
    ***********************/

void print_MR(PM_mapping_results * MR,
              int printPairedLinks);

#ifdef __cplusplus
}
#endif

#endif // PM_BAM_PARSER_H