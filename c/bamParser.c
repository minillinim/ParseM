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
#include "bamParser.h"
#include "pairedLink.h"

// proper linking read is a properly paired, (primary alignment) of the first read in thr pair
#define PM_BAM_FSUPP (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)
#define PM_BAM_FMAPPED (BAM_FMUNMAP | BAM_FUNMAP)

// This function reads a BAM alignment from one BAM file.
int read_bam(void *data,
             bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret = aux->iter? hts_itr_next(aux->fp, aux->iter, b, (hts_readrec_f)(bam_readrec), 0) : bam_read1(aux->fp, b);
    if (!(b->core.flag&BAM_FUNMAP)) {
        if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
        else if (aux->min_len && bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
    }
    return ret;
}


void init_MR(PM_mapping_results * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             int doAverages)
{
    //----
    // initialise the MR object
    //
    int i = 0;
    MR->num_contigs = BAM_header->n_targets;
    MR->num_bams = numBams;
    
    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        // only make these guys if needed
        if(doAverages) {
            MR->total_reads = calloc(MR->num_contigs, sizeof(int*));
            for(i = 0; i < MR->num_contigs; ++i) {
                MR->total_reads[i] = calloc(MR->num_bams, sizeof(int));
            }        
        }
        // always make these guys
        MR->contig_names = calloc(MR->num_contigs, sizeof(char*));
        MR->contig_lengths = calloc(MR->num_contigs, sizeof(float));
        for(i =0; i < MR->num_contigs; ++i) {
            MR->contig_names[i] = strdup(BAM_header->target_name[i]);
            MR->contig_lengths[i] = (float)BAM_header->target_len[i];
        }
    }
    
    cfuhash_table_t *links = cfuhash_new_with_initial_size(30);  
    cfuhash_set_flag(links, CFUHASH_FROZEN_UNTIL_GROWS);
    MR->links = links;
}

void destroy_MR(PM_mapping_results * MR)
{
    //----
    // free the MR object
    //
    int i = 0;
    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        if(MR->total_reads != NULL) {
            for(i = 0; i < MR->num_contigs; ++i) {
                if(MR->total_reads[i] != NULL)
                    free(MR->total_reads[i]);
            }
            free(MR->total_reads);
        }
        if(MR->contig_names != NULL) {
            for(i = 0; i < MR->num_contigs; ++i) {
                if(MR->contig_names[i] != NULL)
                    free(MR->contig_names[i]);
            }
            free(MR->contig_names);
        }
        if(MR->contig_names != NULL)
            free(MR->contig_lengths);
    }
    
    // destroy paired links
    destroyLinks(MR->links);
    cfuhash_clear(MR->links);
    cfuhash_destroy(MR->links);
    
    free(MR);
}

void print_MR(PM_mapping_results * MR, int printPairedLinks) {
    //----
    // print mapping results
    //
    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        if(MR->total_reads != NULL) {
            int i = 0, j = 0;
            for(i = 0; i < MR->num_contigs; ++i) {
                printf("%s\t%0.0f", MR->contig_names[i], MR->contig_lengths[i]);
                for(j = 0; j < MR->num_bams; ++j) {
                    printf("\t%0.3f", (float)MR->total_reads[i][j]/MR->contig_lengths[i]);
                }
                printf("\n");
            }
        }
    }
    if(printPairedLinks) {
        printLinks(MR->links, MR->contig_names);
    }
}

int processBams(int numBams,
                int baseQ,
                int mapQ,
                int minLen, 
                int doLinks,
                int ignoreSuppAlignments,
                int doOutlierCoverage,
                char* bamFiles[],
                PM_mapping_results * MR
) {
    //-----
    // work out coverage depths and also pairwise linkages if asked to do so
    //
    
    // initialize the auxiliary data structures
    const bam_pileup1_t **plp;
    bam_hdr_t *h = 0; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int i = 0, tid = -1, pos, *n_plp;
    int contig_counter = 0; // used for populating mapping results

    // load contig names and BAM index.
    data = calloc(numBams, sizeof(void*)); // data[i] for the i-th input
    int beg = 0, end = 1<<30;  // set the default region
    for (i = 0; i < numBams; ++i) {
        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = bgzf_open(bamFiles[i], "r"); // open BAM
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = minLen;                 // set the qlen filter
        bam_hdr_t *htmp;
        htmp = bam_hdr_read(data[i]->fp);         // read the BAM header ( I think this must be done!)
        if (i == 0) {
            h = htmp; // keep the header of the 1st BAM
        } else bam_hdr_destroy(htmp); // if not the 1st BAM, trash the header

            if (tid >= 0) { // if a region is specified and parsed successfully
                hts_idx_t *idx = bam_index_load(bamFiles[i]);  // load the index
                data[i]->iter = bam_itr_queryi(idx, tid, beg, end); // set the iterator
                hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
            }
    }
    
    // initialise the mapping results struct
    init_MR(MR, h, numBams, 1);
    
    // the core multi-pileup loop
    mplp = bam_mplp_init(numBams, read_bam, (void**)data); // initialization
    n_plp = calloc(numBams, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(numBams, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)
    
    // initialize
    for (i = 0; i < numBams; ++i) {
        MR->total_reads[contig_counter][i] = 0;
    }
    int prev_tid = -1;  // the id of the previous positions tid        
    
    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if(prev_tid != -1 && tid != prev_tid) {  // we've arrived at a new chrom
            for (i = 0; i < numBams; ++i) {
                MR->total_reads[contig_counter+1][i] = 0;  // reset for next chrom
            }
            ++contig_counter;
        }
        
        for (i = 0; i < numBams; ++i) { // base level filters have to go here
            int j, m = 0;
            for (j = 0; j < n_plp[i]; ++j) {
                int bad  = 0;
                // for each read in the pileup
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
                if (p->is_del || p->is_refskip) {++m; ++bad;} // having dels or refskips at tid:pos
                else if (bam_get_qual(p->b)[p->qpos] < baseQ) {++m; ++bad;} // low base quality
                // now we do links if we've been asked to
                if(doLinks && bad == 0) {
                    bam1_core_t core = p->b[0].core;
                    // check to see if this is a proper linking paired read
                    if (p->is_head &&                               // first time we've seen this read
                        (core.flag & BAM_FPAIRED) &&                // read is a paired read
                        (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                        ((core.flag & PM_BAM_FMAPPED) == 0) &&      // both ends are mapped
                        ((core.flag & PM_BAM_FSUPP) == 0) &&        // is primary mapping
                        core.tid != core.mtid) {                    // hits different contigs
                        
                        // looks legit
                        addLink(MR->links,core.tid,                 // contig 1
                                core.mtid,                          // contig 2
                                core.pos,                           // pos 1
                                core.mpos,                          // pos 2
                                ((core.flag&BAM_FREVERSE) != 0),    // 1 == reversed
                                ((core.flag&BAM_FMREVERSE) != 0),   // 0 = agrees
                                i);                                 // bam file ID
                    }
                }
            }
            MR->total_reads[contig_counter][i] += n_plp[i] - m;  // add this positions depth to the total
        }
        prev_tid = tid;
    }
    
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);
    bam_hdr_destroy(h);
    
    for (i = 0; i < numBams; ++i) {
        bgzf_close(data[i]->fp);
        if (data[i]->iter) bam_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data);
    return 0;
}
