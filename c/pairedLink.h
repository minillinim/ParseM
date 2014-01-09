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

#ifndef PM_PAIRED_LINK_H
#define PM_PAIRED_LINK_H

// system includes
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

// cfuhash
#include "cfuhash.h"

/*! @typedef
 @ *a*bstract Structure for storing information about specific links
 @field pos_1 position of read in contig 1
 @field orient_1 orientation of read in contig 1 (1 == reversed)
 @field pos_2 position of read in contig 2
 @field orient_2 orientation of read in contig 2 (1 == reversed)
 @field bam_ID id of the BAM file link originates from
 @field next_link link to next link info struct (linked list)
 */
typedef struct {
    int pos_1;
    int orient_1;
    int pos_2;
    int orient_2;
    int bam_ID;
    struct PM_link_info * next_link;
} PM_link_info;

/*! @typedef
 * @ *a*bstract Structure for storing information about specific links
 * @field cid_1 tid of contig 1 (from BAM header)
 * @field cid_2 tid of contig 2 (cid_1 < cid_2)
 * @field numLinks number of link info structures in linked list
 * @field LI first link info struct for this contig pair
 */
typedef struct {
    int cid_1;
    int cid_2;
    int numLinks;
    PM_link_info * LI;
} PM_link_pair;

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @abstract Add a link info struct to the main link table.
 * 
 * @param  cid_1  tid of contig 1 ( from BAM header )
 * @param  cid_2  tid of contig 2
 * @param  pos_1  position of read in contig 1
 * @param  pos_2  position of read in contig 2
 * @param  orient_1  orientation of read in contig 1
 * @param  orient_2  orientation of read in contig 2
 * @param  bam_ID  id of the BAM file link originates from
 * @return void
 * 
 * @discussion The contigs can be added in any order hovever the pos and orient_1
 * variables must match this ordering. The code will sort cids accordingly.
 */
void addLink(cfuhash_table_t * linkTable,
             int cid_1,
             int cid_2,
             int pos_1,
             int pos_2,
             int orient_1,
             int orient_2,
             int bam_ID);

/*!
 * @abstract Walk along the link info linked list destoying the current node
 * 
 * @param  LI_ptr pointer to the current node
 * @return 1 if should keep walking, 0 is at end of list
 */
int destroyLinkInfo_andNext(PM_link_info** LI_ptr);

/*!
 * @abstract Destroy all links information
 * 
 * @param  linkHash pointer to the hash to be desroyed
 * @return void
 */
void destroyLinks(cfuhash_table_t * linkHash);

/*!
 * @abstract Walk along the link info linked list
 * 
 * @param  LI_ptr pointer to the current node
 * @return 1 if should keep walking, 0 is at end of list
 */
int getNextLinkInfo(PM_link_info** LI_ptr);

        /***********************
        *** PRINTING AND I/O ***
        ***********************/

void printLinks(cfuhash_table_t * linkHash, char ** contigNames);
void printLinkPair(PM_link_pair* LP, char ** contigNames);
void printLinkInfo(PM_link_info* LI);

#ifdef __cplusplus
}
#endif

#endif // PM_PAIRED_LINK_H