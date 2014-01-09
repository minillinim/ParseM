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
#include <stdio.h>
#include <unistd.h>

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"

void makeContigKey(char* keyStore, int cid_1, int cid_2)
{
    if(cid_1 < cid_2) {sprintf(keyStore, "%d,%d", cid_1, cid_2);}
    else {sprintf(keyStore, "%d,%d",cid_2, cid_1);}
}

void addLink(cfuhash_table_t * linkHash,
             int cid_1,
             int cid_2,
             int pos_1,
             int pos_2,
             int orient_1,
             int orient_2,
             int bam_ID
            )
{
    // store the link info
    PM_link_info* LI = (PM_link_info*) calloc(1, sizeof(PM_link_info));
    if(cid_1 < cid_2){
        LI->orient_1 = orient_1;
        LI->orient_2 = orient_2;
        LI->pos_1 = pos_1;
        LI->pos_2 = pos_2;
    }
    else
    {
        LI->orient_1 = orient_2;
        LI->orient_2 = orient_1;
        LI->pos_1 = pos_2;
        LI->pos_2 = pos_1;
    }
    LI->bam_ID = bam_ID;
    
    PM_link_info** next_link_ptr = (PM_link_info**) &LI->next_link;
    
    // see if the key is in the hash already
    char * key = calloc(30, sizeof(char)); // allocate room for the key
    makeContigKey(key, cid_1, cid_2);
    PM_link_pair * base_LP = cfuhash_get(linkHash, key);
    if (base_LP != NULL)
    {
        // exists in the hash -> daisy chain it on
        *next_link_ptr = base_LP->LI;
        base_LP->LI = LI;
        base_LP->numLinks++;
    }
    else 
    {
        // we'll need to build a bit of infrastructure
        // store the contig ids once only
        PM_link_pair * LP = (PM_link_pair*) calloc(1, sizeof(PM_link_pair));
        if(cid_1 < cid_2)
        {
            LP->cid_1 = cid_1;
            LP->cid_2 = cid_2;
        }
        else
        {
            LP->cid_1 = cid_2;
            LP->cid_2 = cid_1;
        }
        LP->LI = LI;    
        LP->numLinks = 1;
        *next_link_ptr = LI; // point to self means end of list
        
        // finally, add the lot to the hash
        cfuhash_put(linkHash, key, LP);
    }
    free(key);
}

int destroyLinkInfo_andNext(PM_link_info** LI_ptr)
{
    PM_link_info* next_link = (PM_link_info*) (*LI_ptr)->next_link;
    if(*LI_ptr ==  next_link) // at the end of the chain
    {
        free(*LI_ptr);
        return 0;
    }
    else // not done yet
    {
        PM_link_info* tmp_link = *LI_ptr;
        free(tmp_link);
        *LI_ptr = next_link;
        return 1;
    }
}

void destroyLinks(cfuhash_table_t * linkHash)
{
    char **keys = NULL;
    size_t *key_sizes = NULL;
    size_t key_count = 0;
    int i = 0;
    keys = (char **)cfuhash_keys_data(linkHash, &key_count, &key_sizes, 0);
    
    for (i = 0; i < (int)key_count; i++) {
        PM_link_pair * base_LP = cfuhash_get(linkHash, keys[i]);
        free(keys[i]);
        PM_link_info* LI = base_LP->LI;
        while(destroyLinkInfo_andNext(&LI));
        free(base_LP);
    }
    free(keys);
    free(key_sizes);    
}

int getNextLinkInfo(PM_link_info** LI_ptr)
{
    PM_link_info* next_link = (PM_link_info*) (*LI_ptr)->next_link;
    if(*LI_ptr ==  next_link) {return 0;}
    else
    {
        *LI_ptr = next_link;
        return 1;
    }
}

void printLinks(cfuhash_table_t * linkHash, char ** contigNames) 
{
    char **keys = NULL;
    size_t *key_sizes = NULL;
    size_t key_count = 0;
    int i = 0;
    keys = (char **)cfuhash_keys_data(linkHash, &key_count, &key_sizes, 0);
    
    for (i = 0; i < (int)key_count; i++) {
        PM_link_pair * LP = cfuhash_get(linkHash, keys[i]);
        free(keys[i]);
        printLinkPair(LP, contigNames);
    }
    free(keys);
    free(key_sizes);    
}

void printLinkPair(PM_link_pair* LP, char ** contigNames) 
{
    printf("===\n(%s, %s, %d links)\n",  contigNames[LP->cid_1], contigNames[LP->cid_2], LP->numLinks);
    PM_link_info* LI = LP->LI;
    do {
        printf("\t");
        printLinkInfo(LI);
        printf("\n");
    } while(getNextLinkInfo(&LI));
}

void printLinkInfo(PM_link_info* LI) 
{
    printf("(%d,%d -> %d,%d, BAM %d)",  LI->pos_1, LI->orient_1, LI->pos_2, LI->orient_2, LI->bam_ID);
}
