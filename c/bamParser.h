//#############################################################################
//
//   __Script__Name__
//   
//   <one line to give the program's name and a brief idea of what it does.>
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

typedef struct {                    // auxiliary data structure
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
    int ** total_reads;             // stores number of mapped reads for contigs
    char ** contig_names;
    float * contig_lengths;
    int num_bams;                   // mainly used so that the free operation works nicely
    int num_contigs;
    cfuhash_table_t * links;        // stores linking information
} PM_mapping_results;

#ifdef __cplusplus
extern "C" {
#endif

int read_bam(void *data,
                bam1_t *b);

void init_MR(PM_mapping_results * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             int doAverages);

void destroy_MR(PM_mapping_results * MR);

void print_MR(PM_mapping_results * MR,
              int printPairedLinks);

int processBams(int numBams,
                int baseQ,
                int mapQ,
                int minLen, 
                int doLinks,
                int ignoreSuppAlignments,
                int doOutlierCoverage,
                char* bamFiles[],
                PM_mapping_results * MR);

#ifdef __cplusplus
}
#endif

#endif // PM_BAM_PARSER_H