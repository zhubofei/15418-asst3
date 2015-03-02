#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "CycleTimer.h"
#include "bfs.h"
#include "graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
#define SPAN 1000

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->alloc_count = count;
    list->present = (int*)malloc(sizeof(int) * list->alloc_count);
    vertex_set_clear(list);
}

void bottom_up_step_omp(
    graph* g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances,
	  int step)
{
    int N = g->num_nodes;
    #pragma omp parallel for
    for (int i=0; i < N; i+=SPAN) {
      int new_outgoings[SPAN];
  		int count = 0;
  		int node_index;

  		for(int j=0; j<SPAN; j++){
  			node_index = i+j;

  			if(node_index < N && distances[node_index] == NOT_VISITED_MARKER){
  				int start_edge = g->incoming_starts[node_index];
  				int end_edge = (node_index < g->num_nodes-1) ? g->incoming_starts[node_index+1] : g->num_edges;

  				for (int p=start_edge; p<end_edge; p++) {
  					int incoming = g->incoming_edges[p];

  					if (distances[incoming] == step-1) {
  						distances[node_index] = step;

  						new_outgoings[count++] = node_index;

  						break;
  					}
  				}
  			}
  		}

  		if(count>0){
        int index;
  			do{
  				index = new_frontier->count;
  			} while(!__sync_bool_compare_and_swap(&new_frontier->count, index, index+count));

  			for(int i=0;i<count;i++){
  				new_frontier->present[index+i] = new_outgoings[i];
  			}
  		}
    }
}

void bfs_bottom_up(graph* graph, solution* sol)
{
  vertex_set list1;
  vertex_set list2;
  vertex_set_init(&list1, graph->num_nodes);
  vertex_set_init(&list2, graph->num_nodes);

  vertex_set* frontier = &list1;
  vertex_set* new_frontier = &list2;

  // initialize all nodes to NOT_VISITED
  // parallel for
  #pragma omp parallel for
  for (int i=0; i<graph->num_nodes; i++) {
      sol->distances[i] = NOT_VISITED_MARKER;
  }

  // setup frontier with the root node
  frontier->present[frontier->count++] = ROOT_NODE_ID;
  sol->distances[ROOT_NODE_ID] = 0;

  int step = 0;

  while (frontier->count != 0) {
      step++;

  #ifdef DEBUG
      double start_time = CycleTimer::currentSeconds();
  #endif

      vertex_set_clear(new_frontier);

      bottom_up_step(graph, frontier, new_frontier, sol->distances, step);

  #ifdef DEBUG
      double end_time = CycleTimer::currentSeconds();
      printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
  #endif

      // swap pointers
      vertex_set* new_outgoings = frontier;
      frontier = new_frontier;
      new_frontier = new_outgoings;
  }

  free(frontier->present);
  free(new_frontier->present);
}


// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    graph* g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances,
    int step)
{
    #pragma omp parallel for
    for (int i=0; i<frontier->count; i++) {

        int node = frontier->present[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes-1) ? g->num_edges : g->outgoing_starts[node+1];

        // attempt to add all neighbors to the new frontier if end_edge > start_edge
        if (end_edge > start_edge) {
            int new_outgoings[end_edge - start_edge];
            int count = 0;

            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                int outgoing = g->outgoing_edges[neighbor];

                if (distances[outgoing] == NOT_VISITED_MARKER) {
                    distances[outgoing] = step;
                    new_outgoings[count++] = outgoing;
                }
            }

            if (count > 0) {
              int index;
              do {
                index = new_frontier->count;
              } while(!__sync_bool_compare_and_swap(&new_frontier->count, index, index+count));

              for(int i=0; i < count; i++)
                new_frontier->present[index + i] = new_outgoings[i];
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(graph* graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    // parallel for
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    frontier->present[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    int step = 0;

    while (frontier->count != 0) {
        step++;

#ifdef DEBUG
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances, step);

#ifdef DEBUG
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* new_outgoings = frontier;
        frontier = new_frontier;
        new_frontier = new_outgoings;
    }

    free(frontier->present);
    free(new_frontier->present);
}
