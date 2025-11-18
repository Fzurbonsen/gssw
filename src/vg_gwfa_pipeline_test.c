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


void test_case1() {

    // set alignment parameters
    int8_t match = 1;
    int8_t mismatch = 1;
    uint8_t gap_open = 1;
    uint8_t gap_extension = 1;
    int8_t full_length_bonus = 0;

    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);



    // set read
    const char* read = "GATCCAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTTTTGTAGTGTCTATAAGTGAACATTTGGCGTGGTTTCAGGCGTAACGTGAAAAAGGAAA";



    // create graph
    gssw_node *node1, *node2;
    node1 = gssw_node_create("Node1", 1, "CATAGAGCAGGTTTGAAACACTCTTTTTGTAGTATCTGGATGTGGACATTTGGATCGCTTTCAGGCCTATGGTGAAAAAGGAAATATCTTCCCATGAAAACTAGACAGAAGCATTCTCAGAAACTTATTTGTGATGTGTGCCCTCAACTGACAGTGTTGAACCTTTGTTTTGATAGAGCAGTTCTGAAACACACTTTTTGTAAAATCTGCAAGAGGATATTTGGATAGCTTTGAGGATTTCGTTGGAAACGGGAATGTCTTCATGTAAACTCTAGACAGAAGCATTCTCAGAAACTGCTTTGGGATGTTTCAATTGAAGTCCCAGTGTTGAACATTCCCTTTCATAGAGCAGGTTTGAAACACTCTTTTTGTACTATCTGGAAGTGGACATTTGGAGCGCTTTCAGGTCTACGGTGAAAAAGGAGATATCTTCCAATAAAAACTAGATAGAAGCAATGTCAGAACTTTTTTCATGATGTATCTACTCAGCAAACAGAGTTGAACCTTTCTTTTGAGAGAGCAGTTTTGAAACACTCCTTTTGTGGAATATGCAAGTGGGTATTAGGCCAGCTTGGAGGATTTCGTTGGAAACGGGAATACGTATAAAAAGCAGACAGCAGCATTGTCAGAAACTACTTTGTGATGTTTGCATTCAAGTCACAGAATTGAACACACCCTTTCACAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTGTAAGTGAACATTTGGATTGATTTCAGGCCTAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTAGACAGAAGCATTCTCAGAAACTTGTTTGTGATGTGTGCCCTCTACTGACAGAGTTGAACCTTTCTTTGCAAAGAGCAGTTTTGAAACACTCTTTTTGTAGAATCTGCAAGAGGATACTTGGATAGCTTTGAGGATTTCTTGGGAAACGGGAATGTCTTCAGATAAACTCTAGACAGAAGCATTCTCAGAAACTTCTTTGGGATGTTTCAATTGAAGTCACAGTGTTGAACATTCCCTTTCA", nt_table, mat);
    node2 = gssw_node_create("Node2", 2, "CAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTATAAGTGAACATTTGGCGTGCTTTCAGGCGTAACGTGAAAAAGGAAATATCTTCCCATAAAAACTAGACAGAAGCATTCTCAGAAACTTGTTCGTGATGTGTGCCCTCTACTGACAGAGTTGAACCTTTCTTTGCAAAGAGCAGCTTTGAAACACACTTTTTGTAGAATCTGCAAGAGGATATTTGGATAGTTTTGAGGATTTCGTTGGAAACGGGTATGTCTTCAGATAAACTCTAGACAGAAGCATTCTCAGAAACTTCTTTGGGATGTTGCATTCAAGTCACAGAGTAGAACATTCCCATTCATAGAGCAGATTTGAAACACTCTTTTTGTAGTATCTGGAAGTGGACATTTGGAGCGCTTTCAGGCCTATGTTGAAAAAGGAAATATCTTCCCATAAAAACTAGACGGAAGCATTCTCAGAAACTTACTTGTGATGTGTTTGCTCAACTAACAGAATTGAACCATCGTTTGGAAGGAGCAGTTTTGAAACACTGTTTTCGTGGAATCTGCAAGTGGATATTTGGCTAGCTTTGAGGATTTCGTTGGAAACGGGATTACATATAAAAAGGAGACAGCAGCATTCTCAGAAACTTCTTTGTGATGTCTGCATTCAAGTCACAGAGTTGAGCATTCCCTTTCATAGAGCAGGTTGGAAACACTCTTTTTGTAGTATCTGGATGAGGACATTTGGAGCGCTTTCAGGCGTATGGTGAAAAAGGAAATATCTTCCCGTAAAAACTAGACAGAAGCATTCTCAGAAATTTATTTGTGATGTGTGCCCTCAACTAACAGAGTTGAACCTTTCTTTTGATAGAGCAGTTTTGAAACACTCTTTTTGTAAAATCTGCAAGAGGATATTTGGATAGCTTGGAGGATTTCATTGCAAACGGGAATGGCTTCATATAAACTCTAGACAGAAGCATTCTCAGAAACTTCGTTGGGATGTTTCGATTGAAGTCCCAGTGTTGAACATTCCCTTTTATAG", nt_table, mat);

    gssw_nodes_add_edge(node1, node2);

    gssw_graph* graph = gssw_graph_create(2);
    gssw_graph_add_node(graph, node1);
    gssw_graph_add_node(graph, node2);



    // perform alignment
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
                                                            full_length_bonus,
                                                            GWFA_CSSWL_INFIX,
                                                            1);

    gssw_print_graph_cigar(&gm1->cigar, stderr);
    fprintf(stderr, "gssw: %i\n", gm1->score);
    gssw_print_graph_cigar(&gm2->cigar, stderr);
    fprintf(stderr, "gwfa: %i\n", gm2->score);

    gssw_graph_mapping_destroy(gm1);
    gssw_graph_mapping_destroy(gm2);
    gssw_graph_destroy(graph);
    free(nt_table);
    free(mat);
}





void test_case2() {

    // set alignment parameters
    int8_t match = 1;
    int8_t mismatch = 1;
    uint8_t gap_open = 1;
    uint8_t gap_extension = 1;
    int8_t full_length_bonus = 0;

    gssw_sse2_disable();

    for (int i = 0; i < 50; ++i) {
        int8_t* nt_table = gssw_create_nt_table();
        int8_t* mat = gssw_create_score_matrix(match, mismatch);



        // set read
        const char* read = "GGCGACAGAGCGAGACTCCGTCTCAAAAAAAAACAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTGGATCATCTCGTATGCCGTCTTCTGCTTGGGCGACAGAGCGAGACT";



        // create graph
        gssw_node *node1, *node2, *node3, *node4, *node5, *node6, *node7, *node8;
        node1 = gssw_node_create("Node1", 1, "TGGTGGCGGGCGCCTGTAGTCCTGGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCGGAGATCGCGCCACTGCACTCCAGCCTGGG", nt_table, mat);
        node2 = gssw_node_create("Node2", 2, "C", nt_table, mat);
        node3 = gssw_node_create("Node3", 3, "A", nt_table, mat);
        node4 = gssw_node_create("Node4", 4, "AACAGAGCGAGACTCCGTCTCAAAAAAAAAAAATTGTTTAAATTGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAATTGTTTAAATTG", nt_table, mat);
        node5 = gssw_node_create("Node5", 5, "G", nt_table, mat);
        node6 = gssw_node_create("Node6", 6, "ACAGAGCGAGACTCCGTCTCAAAAAAAAAAAATTGTTTAAATTGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAATTGTTTAAATTG", nt_table, mat);
        node7 = gssw_node_create("Node7", 7, "A", nt_table, mat);
        node8 = gssw_node_create("Node8", 8, "AAAAAAAAAATGACAGGGTCTTGCTGTCACCCAGGCTGGAGCGTAGTGGCCCTGATCATGATTCACAGTAGCCTCAAACTCCTGGGCTCAAGCAATCTTCCCTCTTCAGCCTCCCAAAGCACTGGAATTACAAGCCTGAGCCACTGCACCTGGCAAGAGGCCATGTTTTTGATCTTGGAATTTCACGGTACCTAAGGCCACATAGTACACCCTCAGGAAAGCAAACAAGTTAATGACAAATTAAAGGTAACCTGTATTTATTGCTGTCGCCTTTCCAATCATGGTGACGTGTCATAGGGTGGGTGCGTTCTTCCTCAGAGACAGACCCCGAGCTGTACCCCAACCTAGTCTTCTCCCAGCTAAATCTCCCTCCCTGGTCACCTCATTCTAGTTCATGAAGTGATTGGCACATTTGCTTAGCTCTACGGCTGCTTCATCACGATTCTTTTAGTCAGTATCTTCACACTAGTATATGAGCTTTCTTA", nt_table, mat);

        // node 1 ->
        gssw_nodes_add_edge(node1, node2);
        gssw_nodes_add_edge(node1, node3);
        // node 2 ->
        gssw_nodes_add_edge(node2, node4);
        gssw_nodes_add_edge(node2, node5);
        // node 3 ->
        gssw_nodes_add_edge(node3, node4);
        gssw_nodes_add_edge(node3, node5);
        // node 4 ->
        gssw_nodes_add_edge(node4, node6);
        // node 5 ->
        gssw_nodes_add_edge(node5, node6);
        // node 6 ->
        gssw_nodes_add_edge(node6, node7);
        gssw_nodes_add_edge(node6, node8);
        // node 7 ->
        gssw_nodes_add_edge(node7, node8);

        gssw_graph* graph = gssw_graph_create(8);
        gssw_graph_add_node(graph, node1);
        gssw_graph_add_node(graph, node2);
        gssw_graph_add_node(graph, node3);
        gssw_graph_add_node(graph, node4);
        gssw_graph_add_node(graph, node5);
        gssw_graph_add_node(graph, node6);
        gssw_graph_add_node(graph, node7);
        gssw_graph_add_node(graph, node8);


        // // perform alignment
        // gssw_graph_fill_pinned(graph,
        //                         read,
        //                         nt_table,
        //                         mat,
        //                         gap_open,
        //                         gap_extension,
        //                         full_length_bonus,
        //                         full_length_bonus,
        //                         15,
        //                         2,
        //                         1);

        // gssw_graph_mapping* gm1 = gssw_graph_trace_back(graph,
        //                                                 read,
        //                                                 strlen(read),
        //                                                 nt_table,
        //                                                 mat,
        //                                                 gap_open,
        //                                                 gap_extension,
        //                                                 full_length_bonus,
        //                                                 full_length_bonus);

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
                                                                full_length_bonus,
                                                                GWFA_CSSWL_INFIX,
                                                                0);

        if (i == 1 || i == 0) {
            gssw_print_graph_cigar(&gm2->cigar, stderr);
            fprintf(stderr, "gwfa: %i\n", gm2->score);
        }

        // gssw_graph_mapping_destroy(gm1);
        gssw_graph_mapping_destroy(gm2);
        gssw_graph_destroy(graph);
        free(nt_table);
        free(mat);
    }
}



void test_case3() {

    // set alignment parameters
    int8_t match = 1;
    int8_t mismatch = 1;
    uint8_t gap_open = 1;
    uint8_t gap_extension = 1;
    int8_t full_length_bonus = 0;

    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);



    // set read
    const char* read = "TTCAAGCAGAAGACGGCATACGAGATGATCCAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTTTGTGATGTGTGCATTCAACTCACAGAGTTCAA";



    // create graph
    gssw_node *node1;
    node1 = gssw_node_create("Node1", 1, "GCTTTGATGCCTATGGTGGAAAAGGAAATATCCGCCCATAAAAACTAGACAGCAGCATTCTCAGAAAGTTGTTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACCTTTCCTTTGATTGAGCAGTTTTGAAAAAGTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCAGTTTGAGGCCTATGGTGTAAAAGGAAATATCTTCACATGAAAACTAGACAGAAGCATTCTCAGAAACTTCTTTGTGATGAGTTCATTCAATTCACATAGTTGAACATTTCTTTTGATAGAGTAGTTTTGAAACACTCTTTCTGTAGAATCTACAAGTGGATATTTGGAGCACATTGAAGCCTATGATGGAAAAGGAAATATCTTCACATACAAACTAGACAGAAGCATTCTCAGAAACTTCTTTGTGATAAGTGCATTCAACTCACAGAGTCGAACCTTTCTGTTGATAGAGCAGTTTTAAATCACTCTTTTTCTAGAATCTGAAAGTGGATATTTGGAGAGCTTTGAGGCCTATGGTGGAAAAGGAAATACCTACGCATAAAAACTATGCGGAAGCATTCTCAGAAATATCTTTGTGATGAGTGCATTCAACTCACAGAGTTGAACATTTATGTTGATAGAGGAGTTTTAAAACACTCTTTTTCAGGAATCTGAAAGTGGATATTTGGAGCGCTTTGAGGCCTATGGTGGAAAAGGAAACACCTTCACAAAAAAAACTAGAGCAGAATCATTCTCAGGAACTTCTTTGTGATGTGTGCATTCAACTCACAGAGTTGAACCTTTTTATTTGATAGAGCAGTTTTGAAACACTATTTTTGTACAATCTGCGGTTGGATATTTGGAGCGCTTTGATGCCTATGGTGGAAAACGAAATATCCGCACATAAAATCTAGACAGCAGCATTCTCAGAAACTTGTTAGTGTTGTGTGCATTCAGCTCACAGAGTTGAACCTTTCCTTTGATTGAGCAGTTTTGAAATAGTCTTTTTGTAGAATCCACAAGTGGATATTT", nt_table, mat);



    gssw_graph* graph = gssw_graph_create(1);
    gssw_graph_add_node(graph, node1);



    // perform alignment
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

    gssw_print_graph_cigar(&gm1->cigar, stderr);
    fprintf(stderr, "gssw: %i\n", gm1->score);

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
                                                            full_length_bonus,
                                                            GWFA_CSSWL_INFIX,
                                                            1);

    gssw_print_graph_cigar(&gm2->cigar, stderr);
    fprintf(stderr, "gwfa: %i\n", gm2->score);

    gssw_graph_mapping_destroy(gm1);
    gssw_graph_mapping_destroy(gm2);
    gssw_graph_destroy(graph);
    free(nt_table);
    free(mat);
}



void test_case4() {

    // set alignment parameters
    int8_t match = 1;
    int8_t mismatch = 1;
    uint8_t gap_open = 1;
    uint8_t gap_extension = 1;
    int8_t full_length_bonus = 0;

    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);


                    
    // set read      // TTCAAGCAGAAGACGATCAAGTGACTGGAGTTCAGACGTGTGTGTAGTGC
    const char* read = "TTCAAGCAGAAGACGATCAAGTGACTGGAGTTCAGACGTGTGTGTAGTGC";



    // create graph
    gssw_node *node1;                  // TTCAAGCAGAAGACTCACAGCAGGTTTGTAGTGC
    node1 = gssw_node_create("Node1", 1, "TTCAAGCAGAAGACTCACAGCAGGTTTGTAGTGC", nt_table, mat);



    gssw_graph* graph = gssw_graph_create(1);
    gssw_graph_add_node(graph, node1);



    // perform alignment
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

    // gssw_graph_print_score_matrices(graph, read, strlen(read), stdout);

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
                                                            full_length_bonus,
                                                            GWFA_CSSWL_INFIX,
                                                            1);

    gssw_print_graph_cigar(&gm1->cigar, stderr);
    fprintf(stderr, "offset: %i\n", gm1->position);
    fprintf(stderr, "gssw: %i\n", gm1->score);
    gssw_print_graph_cigar(&gm2->cigar, stderr);
    fprintf(stderr, "offset: %i\n", gm2->position);
    fprintf(stderr, "gwfa: %i\n", gm2->score);

    gssw_graph_mapping_destroy(gm1);
    gssw_graph_mapping_destroy(gm2);
    gssw_graph_destroy(graph);
    free(nt_table);
    free(mat);
}



void test_case5() {

    // set alignment parameters
    int8_t match = 1;
    int8_t mismatch = 1;
    uint8_t gap_open = 1;
    uint8_t gap_extension = 1;
    int8_t full_length_bonus = 0;

    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);



    // set read
    const char* read = "CCCCTTTCCCC";



    // create graph
    gssw_node *node1;
    node1 = gssw_node_create("Node1", 1, "CCCCGGGCCCC", nt_table, mat);



    gssw_graph* graph = gssw_graph_create(1);
    gssw_graph_add_node(graph, node1);



    // perform alignment
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

    gssw_print_graph_cigar(&gm1->cigar, stderr);
    fprintf(stderr, "offset: %i\n", gm1->position);
    fprintf(stderr, "gssw: %i\n", gm1->score);

    // gssw_graph_print_score_matrices(graph, read, strlen(read), stdout);

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
                                                            full_length_bonus,
                                                            GWFA_CSSWL_INFIX,
                                                            1);

    gssw_print_graph_cigar(&gm2->cigar, stderr);
    fprintf(stderr, "offset: %i\n", gm2->position);
    fprintf(stderr, "gwfa: %i\n", gm2->score);

    gssw_graph_mapping_destroy(gm1);
    gssw_graph_mapping_destroy(gm2);
    gssw_graph_destroy(graph);
    free(nt_table);
    free(mat);
}



int main() {
    fprintf(stderr, "======================================\n");
    fprintf(stderr, "====Running vg_gwfa_pipeline_tests====\n");
    fprintf(stderr, "======================================\n");


    test_case3();

    fprintf(stderr, "=============Run complete!============\n");
    return 0;
}