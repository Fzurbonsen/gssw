/*

    gssw:
    vg_hgwfa_pipeline_wrapper.h
    Wrapper for cpp.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include "gssw.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GWFA_EDLIB_PREFIX,
    GWFA_EDLIB_INFIX,
    GWFA_CSSWL_PREFIX,
    GWFA_CSSWL_INFIX
} vg_gwfa_pipeline_algorithm_type_e;

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
                                                    int8_t end_full_length_bonus,
                                                    vg_gwfa_pipeline_algorithm_type_e algorithm_type,
                                                    int8_t print_debug);

void test();

#ifdef __cplusplus
}
#endif

