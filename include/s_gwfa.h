/*

    s_gwfa.h
    Header file that holds the definitons for s_gwfa.c
    Author: Frederic zur Bonsen <fzurbonsen@ethz.ch>

*/

#ifndef S_GWFA_H
#define S_GWFA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <sys/time.h>


/*

    Graphs need to have ID's corresponding to their position in the node list.

*/


// Structs:

// Forward declarations
struct s_gwfa_node_t;
struct s_gwfa_edge_t;
struct s_gwfa_graph_t;
struct s_gwfa_path_t;

// Node struct
typedef struct s_gwfa_node_t {
    int32_t id; // Id of the node
    int32_t idx; // Index of the node
    int32_t size; // Size of the node
    char* seq; // Sequence stored in the node
    uint16_t n_edges; // Number of outgoing edges
    struct s_gwfa_edge_t** edges; // Pointer to list of outgoing edges
} s_gwfa_node_t;

// Edge struct
typedef struct s_gwfa_edge_t {
    s_gwfa_node_t* start; // Start node of the edge
    s_gwfa_node_t* end; // End node of the edge
    int32_t weight; // Weight of the edge
} s_gwfa_edge_t;

// Graph struct
typedef struct s_gwfa_graph_t {
    uint16_t size;
    s_gwfa_node_t** nodes;
    s_gwfa_node_t* top_node;
    struct s_gwfa_path_t* top_nodes;
} s_gwfa_graph_t;

typedef struct s_gwfa_path_t {
    uint32_t size;
    uint32_t alloc_size;
    s_gwfa_node_t** nodes;

    // for easy deletion
    struct s_gwfa_path_t* next;
} s_gwfa_path_t;


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

// Functions:
// Nodes
s_gwfa_node_t* s_gwfa_node_create(int32_t id,
                                    int32_t size,
                                    const char* seq);

// Edges
s_gwfa_edge_t* s_gwfa_edge_create(s_gwfa_node_t* start,
                                    s_gwfa_node_t* end,
                                    int32_t weight);

// Graphs
s_gwfa_graph_t* s_gwfa_graph_create();
void s_gwfa_graph_add_node(s_gwfa_graph_t* graph, 
                            s_gwfa_node_t* node);
s_gwfa_node_t* s_gwfa_graph_find_top_node(s_gwfa_graph_t* graph);
s_gwfa_path_t* s_gwfa_graph_find_top_nodes(s_gwfa_graph_t* graph);

// Path
s_gwfa_path_t* s_gwfa_path_create();
void s_gwfa_path_push(s_gwfa_path_t* path, s_gwfa_node_t* node);
s_gwfa_path_t* s_gwfa_path_copy(s_gwfa_path_t* path);


// Clean-up
void s_gwfa_edge_destroy(s_gwfa_edge_t* edge);
void s_gwfa_node_destroy(s_gwfa_node_t* node);
void s_gwfa_graph_destroy(s_gwfa_graph_t* graph);
void s_gwfa_path_destroy(s_gwfa_path_t* path);


// ed_wfa

typedef struct wfa_diag_t {
    int32_t idx; // index of the diagonal
    int32_t v_idx; // Distance of the furthest vertex on the diagonal
    int32_t ov_idx; // Old v_idx to keep track during wfa_expand
    int32_t wf; // Wave front
    int8_t done; // Indicator as to wether a diagonal has been updated already
    struct wfa_diag_t* next; // Pointer for queue
} wfa_diag_t;

typedef struct wfa_diag_queue_t {
    size_t size; // queue size
    wfa_diag_t* tail; // Pointer to the tail of the queue
    wfa_diag_t* head; // Pointer to the head of the queue
} wfa_diag_queue_t;

typedef struct wfa_diag_list_t {
    size_t n_u_diags; // Number of upper diags
    size_t n_l_diags; // Number of lower diags
    size_t alloc_u; // Allocated memory for upper diags
    size_t alloc_l; // Allocated memory for lower diags
    wfa_diag_t* z_diag; // Index zero diagonal
    wfa_diag_t** u_diags; // Diags with pos index
    wfa_diag_t** l_diags; // Diags with neg index
} wfa_diag_list_t;

int32_t wfa_ed(const char* seq, int32_t seq_len, const char* ref, int32_t ref_len);
int32_t wfa_ed_infix(const char* seq, int32_t seq_len, const char* ref, int32_t ref_len);
int32_t wfa_ed_prefix(const char* seq, int32_t seq_len, const char* ref, int32_t ref_len);


// g_wfa

typedef struct g_wfa_diag_t {
    int32_t idx; // index of the diagonal
    int32_t v_idx; // Distance of the furthest vertex on the diagonal
    int32_t ov_idx; // Old v_idx to keep track during wfa_expand
    int32_t n_idx; // node index
    int8_t* is_q;
    struct g_wfa_diag_t** next; // Pointer for queue
    struct s_gwfa_path_t* path; // Pointer to path
} g_wfa_diag_t;

typedef struct g_wfa_diag_queue_t {
    size_t size; // queue size
    g_wfa_diag_t* tail; // Pointer to the tail of the queue
    g_wfa_diag_t* head; // Pointer to the head of the queue
    int8_t q_idx;
} g_wfa_diag_queue_t;

typedef struct g_wfa_diag_list_t {
    int32_t n_idx; // Index of the node that the list corresponds to
    size_t n_u_diags; // Number of upper diags
    size_t n_l_diags; // Number of lower diags
    size_t alloc_u; // Allocated memory for upper diags
    size_t alloc_l; // Allocated memory for lower diags
    g_wfa_diag_t* z_diag; // Index zero diagonal
    g_wfa_diag_t** u_diags; // Diags with pos index
    g_wfa_diag_t** l_diags; // Diags with neg index
    s_gwfa_path_t* path_head; // Pointer to linked list of all the paths
} g_wfa_diag_list_t;

typedef struct g_wfa_diag_list_set_t {
    size_t n; // Number of nodes, i.e. number of diag lists
    g_wfa_diag_list_t** dls; // List of diag lists, ordered by node index
} g_wfa_diag_list_set_t;

int32_t g_wfa_ed(const char* seq, int32_t seq_len, s_gwfa_graph_t* graph, s_gwfa_path_t** final_path);
int32_t g_wfa_ed_infix(const char* seq, int32_t seq_len, s_gwfa_graph_t* graph, s_gwfa_path_t** final_path);

#ifdef __cplusplus
}
#endif	// __cplusplus

#endif // S_GWFA_H