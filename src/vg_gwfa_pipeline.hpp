/*

    projectA:
    vg_hgwfa_pipeline.hpp
    This file holds the definitions for the pipeline from vg to gwfa.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

#include "gssw.h"

#include "gwfa/gwfa.h"
#include "gwfa/gfa.h"
#include "gwfa/gfa-priv.h"
#include "gwfa/ketopt.h"
#include "gwfa/kalloc.h"
#include "gwfa/kseq.h"

#include "edlib/edlib.h"

#include "csswl/ssw.h"


using namespace std;


#ifndef PROJECTA_VG_GWFA_PIPELINE_HPP
#define PROJECTA_VG_GWFA_PIPELINE_HPP

#ifdef __cplusplus
extern "C" {
#endif
    void notice_me();
#ifdef __cplusplus
}
#endif

// Class to handle the vg (gssw) alignment process with gwfa.
class ProjectA_VG_GWFA_Aligner {
    
    private:
        // data from vg
        gssw_graph* vg_graph;
        const char* read;
        int8_t* nt_table;
        int8_t* mat;
        uint8_t gap_open;
        uint8_t gap_extension;

        // data for gwfa
        void* km;
        gwf_graph_t* gwfa_graph;
        int32_t ql;
        // const char* q; already covered by read in vg
        int32_t v0;
        int32_t v1;
        uint32_t max_lag;
        int32_t traceback;

        // map to map gwfa nodes to vg nodes
        unordered_map<gssw_node*, int32_t> node_map1;
        unordered_map<int32_t, gssw_node*> node_map2;

        // internal values
        int32_t score;
        gwf_path_t path;
        string reference;
        string cigar;

        // metadata about process for print function
        bool done_graph;
        bool done_align_s2g;
        bool done_path_to_seq;
        bool done_align_s2s;
        bool done_all;

        // output
        gssw_graph_mapping* gm;


    // constructor/destructor
    public:
        ProjectA_VG_GWFA_Aligner(gssw_graph* vg_grpah,
                                const char* read,
                                int8_t* nt_table,
                                int8_t* mat,
                                uint8_t gap_open,
                                uint8_t gap_extension);
        ~ProjectA_VG_GWFA_Aligner();
    
    // methods
    public:
        void align_edlib(int32_t do_traceback);
        void align_edlib_infix(int32_t do_traceback);
        void align_csswl(int32_t do_traceback);
        void align_csswl_infix(int32_t do_traceback);
        void print(FILE* file);
        gssw_graph_mapping* graph_mapping();

    private:
        void _align_ed(); // align with gwfa edit distance, base algorithm (prefix mode)
        void _align_ed_infix(); // align with gwfa edit distance in infix mode
        void _path_to_seq(); // method to build the reference seqeuence from the path
        void _align_edlib_global(); // align with the edlib s2s algorithm in the global mode
        void _align_edlib_prefix(); // align with the edlib s2s algorithm in the prefix mode
        void _align_edlib_infix(); // align with the edlib s2s algorithm in the infix mode
        void _align_csswl(); // align with the csswl s2s algorithm
        void _print_graph(FILE* file); // method to print the graphs
        void _print_path(FILE* file); // method to print the gwfa path
        void _print_graph_cigar(FILE* file); // method to print the graph-CIGAR
        void _gssw_to_gwfa(); // transforms the gssw input data into the gwfa data structures
        void _cigar_to_gssw(); // parses a CIGAR and stores it into the gssw data struct

        
};

#endif // PROJECTA_VG_GWFA_PIPELINE_HPP