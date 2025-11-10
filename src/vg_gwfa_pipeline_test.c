/*
 * tests for the vg_gwfa_pipeline
*/

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "gssw.h"

#include "vg_gwfa_pipeline_wrapper.h"


int main() {
    fprintf(stderr, "======================================\n");
    fprintf(stderr, "====Running vg_gwfa_pipeline_tests====\n");
    fprintf(stderr, "======================================\n");


    int8_t match = 1;
    int8_t mismatch = 4;
    uint8_t gap_open = 6;
    uint8_t gap_extension = 1;
    int8_t full_length_bonus = 0;

    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);

    const char* read = "AAAAA";

    gssw_node* node;
    node = gssw_node_create("Node", 1, "AAAAAA", nt_table, mat);

    gssw_graph* graph = gssw_graph_create(1);
    gssw_graph_add_node(graph, node);

    gssw_graph_fill_pinned(graph,
                            read,
                            nt_table,
                            mat,
                            gap_open,
                            gap_extension,
                            full_length_bonus,
                            full_length_bonus,
                            15,
                            2,
                            1);

    gssw_graph_mapping* gm1 = gssw_graph_trace_back(graph,
                                                    read,
                                                    strlen(read),
                                                    nt_table,
                                                    mat,
                                                    gap_open,
                                                    gap_extension,
                                                    full_length_bonus,
                                                    full_length_bonus);

    gssw_graph_mapping* gm2 = gwfa_graph_align_trace_back(graph,
                                                            0,
                                                            0,
                                                            false,
                                                            read,
                                                            NULL,
                                                            strlen(read),
                                                            NULL,
                                                            0,
                                                            nt_table,
                                                            mat,
                                                            gap_open,
                                                            gap_extension,
                                                            full_length_bonus,
                                                            full_length_bonus);

    gssw_print_graph_cigar(&gm1->cigar, stderr);
    fprintf(stderr, "gssw: %i\n", gm1->score);
    gssw_print_graph_cigar(&gm2->cigar, stderr);
    fprintf(stderr, "gwfa: %i\n", gm2->score);
}