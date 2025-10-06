/*

    gssw:
    vg_hgwfa_pipeline.cpp
    This file holds the implementation for the pipeline from vg to gwfa.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <cstring>
#include <cstdint>

#include "vg_gwfa_pipeline.hpp"
#include "vg_gwfa_pipeline_wrapper.h"


// constructor
ProjectA_VG_GWFA_Aligner::ProjectA_VG_GWFA_Aligner(gssw_graph* vg_graph,
                                                    const char* read,
                                                    int8_t* nt_table,
                                                    int8_t* mat,
                                                    uint8_t gap_open,
                                                    uint8_t gap_extension)
    :   vg_graph(vg_graph),
        read(read),
        nt_table(nt_table),
        mat(mat),
        gap_open(gap_open),
        gap_extension(gap_extension),
        v0(0), // defines the first node as the start node for the alignment
        v1(-1), // no end node
        max_lag(0), // no max lag
        ql(strlen(read)),
        km(::km_init()),
        done_graph(false),
        done_align_s2g(false),
        done_path_to_seq(false),
        done_align_s2s(false),
        done_all(false),
        gm(::gssw_graph_mapping_create())
{
    _gssw_to_gwfa();
}


// destructor
ProjectA_VG_GWFA_Aligner::~ProjectA_VG_GWFA_Aligner() {
    // free gwfa_graph
    if (gwfa_graph) {
        gwf_cleanup(km, gwfa_graph);
        int32_t i;
        free(gwfa_graph->len); free(gwfa_graph->seq); free(gwfa_graph->arc); free(gwfa_graph->src); free(gwfa_graph);
        gwfa_graph = nullptr;
    }
    km_destroy(km);
}


// print the gwfa graph
void ProjectA_VG_GWFA_Aligner::_print_graph(FILE* file) {
    for (int i = 0; i < gwfa_graph->n_vtx; ++i) {
        fprintf(file, "S\t%i\t%s\n", node_map2[i]->id, node_map2[i]->seq);
    }
    for (int i = 0; i < gwfa_graph->n_arc; ++i) {
        fprintf(file, "L\t%i\t+\t%i\t+\t0M\n", node_map2[gwfa_graph->arc[i].a >> 32]->id, node_map2[gwfa_graph->arc[i].a]->id);
    }
}

// print the gwfa path
void ProjectA_VG_GWFA_Aligner::_print_path(FILE* file) {
    for (int i = 0; i < path.nv; ++i) {
        fprintf(file, "[%i]->", node_map2[path.v[i]]->id);
    }
    fprintf(file, "\n");
    fprintf(file, "path score: %i\n", path.s);
    fprintf(file, "end node: %i\n", path.end_v);
    fprintf(file, "end offset: %i\n", path.end_off);
}

// print the gwfa path
void ProjectA_VG_GWFA_Aligner::_print_graph_cigar(FILE* file) {

}


// method to transform gssw data structures to gwfa data structures
void ProjectA_VG_GWFA_Aligner::_gssw_to_gwfa() {

    // initialize gwfa graph
    GFA_CALLOC(gwfa_graph, 1);
    gwfa_graph->n_vtx = vg_graph->size;
    
    GFA_MALLOC(gwfa_graph->len, gwfa_graph->n_vtx);
    GFA_MALLOC(gwfa_graph->src, gwfa_graph->n_vtx);
    GFA_MALLOC(gwfa_graph->seq, gwfa_graph->n_vtx);

    // add nodes
    int32_t n_arc = 0;
    for (int i = 0; i < vg_graph->size; ++i) {
        gssw_node* node = vg_graph->nodes[i];
        gwfa_graph->seq[i] = node->seq;
        gwfa_graph->len[i] = node->len;
        n_arc += node->count_prev; // count the number of edges
        node_map1[node] = i;
        node_map2[i] = node;
    }

    gwfa_graph->n_arc = n_arc;
    GFA_MALLOC(gwfa_graph->arc, gwfa_graph->n_arc);

    // add edges
    for (int i = 0, k = 0; i < vg_graph->size; ++i) {
        gssw_node* node = vg_graph->nodes[i];
        for (int j = 0; j < node->count_prev; ++j, ++k) {
            uint32_t origin = node_map1[node->prev[j]];
            gwfa_graph->arc[k].a = (uint64_t)origin<<32 | i;
            gwfa_graph->arc[k].o = 0;
        }
    }
    done_graph = true;
}


// method to build a sequence from the path
void ProjectA_VG_GWFA_Aligner::_path_to_seq() {
    for (int i = 0; i < path.nv; ++i) {
        reference += gwfa_graph->seq[path.v[i]];
    }
    done_path_to_seq = true;
}


// method to transform the CIGAR string into the gssw graph-CIGAR
void ProjectA_VG_GWFA_Aligner::_cigar_to_gssw() {

    // flatten the CIGAR to make it easier to handle
    string f_cigar; // flattened CIGAR
    int num = 0;

    for (const char* p = cigar.c_str(); *p; ++p) {
        if (isdigit(*p)) {
            num = num * 10 + (*p - '0');
        } else {
            if (num == 0) num = 1;
            f_cigar.append(num, *p);
            num = 0;
        }
    }

    // iterate over the gssw nodes to add into the
    int32_t ref_pos = gm->position; // offset in the first node
    int32_t cigar_idx = 0;

    // create graph CIGAR struct for gssw
    gm->cigar.length = path.nv;
    gm->cigar.elements = (gssw_node_cigar*)malloc(path.nv * sizeof(gssw_node_cigar));
    int32_t local_score = 0;

    int match = mat[0];
    int mismatch = mat[0];
    int insertion = gap_open;
    int deletion = gap_open;

    int counter = ql; // counter to ensure that all of the read is aligned

    // iterate over all the nodes in the path to assign the corresponding cigar
    for (int i = 0; i < path.nv; ++i) {
        gssw_node* node = node_map2[path.v[i]]; // find the node with the help of the node map
        gssw_node_cigar nc;
        nc.node = node;

        gssw_cigar* g_cigar = (gssw_cigar*)calloc(1, sizeof(gssw_cigar));
        int32_t node_size = node->len - ref_pos;

        // we go through the node and assign the CIGAR elements
        while (node_size) {
            if (cigar_idx >= f_cigar.size()
                || !(f_cigar[cigar_idx] == 'M' || f_cigar[cigar_idx] == 'I' || f_cigar[cigar_idx] == 'D' || f_cigar[cigar_idx] == '=' || f_cigar[cigar_idx] == 'X')) {
                break;
            }
            if (f_cigar[cigar_idx] == 'M') {
                node_size--;
                local_score += match;
                counter--;
                gssw_cigar_push_back(g_cigar, f_cigar[cigar_idx], 1);
            } else if (f_cigar[cigar_idx] == 'D') {
                node_size--;
                local_score += deletion;
                gssw_cigar_push_back(g_cigar, f_cigar[cigar_idx], 1);
            } else if (f_cigar[cigar_idx] == 'I') {
                local_score += insertion;
                counter--;
                gssw_cigar_push_back(g_cigar, f_cigar[cigar_idx], 1);
            } else if (f_cigar[cigar_idx] == '=') {
                node_size--;
                local_score += match;
                counter--;
                gssw_cigar_push_back(g_cigar, 'M', 1);
            } else if (f_cigar[cigar_idx] == 'X') {
                node_size--;
                local_score += mismatch;
                counter--;
                gssw_cigar_push_back(g_cigar, 'M', 1);
            }
            cigar_idx++;
        }

        // if we are in the last cycle: ensure that all of the read is aligned
        if (i+1 == path.nv) {
            for (; counter > 0; --counter) {
                gssw_cigar_push_back(g_cigar, 'I', 1);
            }
        }

        nc.cigar = g_cigar;
        ref_pos = 0;
        gm->cigar.elements[i] = nc;
        gm->score = local_score;
    }

    done_all = true;
}


// method to perform edlib global alignment
void ProjectA_VG_GWFA_Aligner::_align_edlib_global() {
    EdlibAlignResult result = edlibAlign(read,
                                        ql,
                                        reference.c_str(),
                                        reference.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    char* edlib_cigar = edlibAlignmentToCigar(result.alignment,
                                        result.alignmentLength,
                                        EDLIB_CIGAR_STANDARD);
    cigar = edlib_cigar;
    free(edlib_cigar);
    gm->score = -result.editDistance;
    gm->position = result.startLocations[0];
    edlibFreeAlignResult(result);
    done_align_s2s = true;
}


// method to perform edlib prefix alignment
void ProjectA_VG_GWFA_Aligner::_align_edlib_prefix() {
    EdlibAlignResult result = edlibAlign(read,
                                        ql,
                                        reference.c_str(),
                                        reference.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    char* edlib_cigar = edlibAlignmentToCigar(result.alignment,
                                        result.alignmentLength,
                                        EDLIB_CIGAR_STANDARD);
    cigar = edlib_cigar;
    free(edlib_cigar);
    gm->score = -result.editDistance;
    gm->position = result.startLocations[0];
    edlibFreeAlignResult(result);
    done_align_s2s = true;
}


// method to perform edlib infix alignment
void ProjectA_VG_GWFA_Aligner::_align_edlib_infix() {
    EdlibAlignResult result = edlibAlign(read,
                                        ql,
                                        reference.c_str(),
                                        reference.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    char* edlib_cigar = edlibAlignmentToCigar(result.alignment,
                                        result.alignmentLength,
                                        EDLIB_CIGAR_STANDARD);
    cigar = edlib_cigar;
    free(edlib_cigar);
    gm->score = -result.editDistance;
    gm->position = result.startLocations[0];
    edlibFreeAlignResult(result);
    done_align_s2s = true;
}


// method to align with ed_infix
void ProjectA_VG_GWFA_Aligner::_align_ed() {

    // index the graph
    ::gwf_ed_index(km, gwfa_graph);

    // perform the alignment
    score = ::gwf_ed(km,
                    gwfa_graph,
                    ql,
                    read,
                    v0,
                    v1,
                    max_lag,
                    traceback,
                    &path);
    done_align_s2g = true;
}


// public method to align with the edlib algorithm in prefix mode
void ProjectA_VG_GWFA_Aligner::align_edlib(int32_t do_traceback) {
    
    if (!(do_traceback == 0 || do_traceback == 1 || do_traceback == 2)) {
        cerr << "[projectA::vg_to_gwfa_pipeline]error: invalid traceback mode!" << endl;
        cerr << "\t" << do_traceback << " is not an allowed traceback mode. Choose one of the following:" << endl
                                                                << "\t0: perform no traceback" << endl
                                                                << "\t1: perform granular traceback" << endl
                                                                << "\t2: perform full traceback in gwfa" << endl;
        exit(1);
    }

    traceback = do_traceback;
    _align_ed();
    _path_to_seq();
    _align_edlib_prefix();
    _cigar_to_gssw();
}


// public method to align with the edlib algorithm in infix mode
void ProjectA_VG_GWFA_Aligner::align_edlib_infix(int32_t do_traceback) {

    if (!(do_traceback == 0 || do_traceback == 1 || do_traceback == 2)) {
        cerr << "[projectA::vg_to_gwfa_pipeline]error: invalid traceback mode!" << endl;
        cerr << "\t" << do_traceback << " is not an allowed traceback mode. Choose one of the following:" << endl
                                                                << "\t0: perform no traceback" << endl
                                                                << "\t1: perform granular traceback" << endl
                                                                << "\t2: perform full traceback in gwfa" << endl;
        exit(1);
    }

    traceback = do_traceback;
    _align_ed();
    _path_to_seq();
    _align_edlib_infix();
    _cigar_to_gssw();
}


void ProjectA_VG_GWFA_Aligner::align_csswl(int32_t do_traceback) {} // ToDo
void ProjectA_VG_GWFA_Aligner::align_csswl_infix(int32_t do_traceback) {} // ToDo


// public method to print the contents of the class
void ProjectA_VG_GWFA_Aligner::print(FILE* file) {
    fprintf(file, "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
    fprintf(file, "Printing the class content:\n");
    if (done_graph) {
        fprintf(file, "\ngraph(1/1):\n");
        _print_graph(file);
    } else {
        fprintf(file, "\ngraph(0/1):\n");
    }

    if (done_align_s2g) {
        fprintf(file, "\npath(1/1):\n");
        _print_path(file);
    } else {
        fprintf(file, "\npath(0/1):\n");
    }

    if (done_path_to_seq) {
        fprintf(file, "\nreference(1/1):\n");
        fprintf(file, "%s\n", reference.c_str());
    } else {
        fprintf(file, "\nreference(0/1):\n");
    }

    if (done_align_s2s) {
        fprintf(file, "\nCIGAR(1/1):\n");
        fprintf(file, "%s\n", cigar.c_str());
    } else {
        fprintf(file, "\nCIGAR(0/1):\n");
    }

    if (done_all) {
        fprintf(file, "\ngraph-CIGAR(1/1):\n");
        gssw_print_graph_cigar(&gm->cigar, file);
    } else {
        fprintf(file, "\ngraph-CIGAR(0/1):\n");
    }
    fprintf(file, "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
}


// method to retrieve gssw graph mapping
gssw_graph_mapping* ProjectA_VG_GWFA_Aligner::graph_mapping() {
    return gm;
}


#ifdef __cplusplus
extern "C" {
#endif

// wrapper function to allow gssw.c to call the traceback
gssw_graph_mapping* gwfa_graph_align_trace_back(gssw_graph* graph,
                                                    int32_t doing_pinning,
                                                    int32_t num_tracebacks,
                                                    int32_t find_internal_node_alts,
                                                    const char* read,
                                                    const char* qual,
                                                    int32_t readLen,
                                                    gssw_node** pinning_nodes,
                                                    int32_t num_pinning_nodes,
                                                    int8_t* nt_table,
                                                    int8_t* score_matrix,
                                                    uint8_t gap_open,
                                                    uint8_t gap_extension,
                                                    int8_t start_full_length_bonus,
                                                    int8_t end_full_length_bonus) {

    // work with class
    ProjectA_VG_GWFA_Aligner aligner(graph,
                                    read,
                                    nt_table,
                                    score_matrix,
                                    gap_open,
                                    gap_extension);

    aligner.align_edlib(1);
    return aligner.graph_mapping();
}

void test() {
    cerr << "error" << endl;
    exit(1);
}

#ifdef __cplusplus
}
#endif