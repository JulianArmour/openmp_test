#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "sys\timeb.h"

typedef struct Node {
    // the id of this node (used for total ordering of nodes)
    int id;
    // color of the node
    int color;
    // the number of neighbors
    int nNeighbors;
    // contains the ids of the neighbors
    int *neighbors;
} Node;


/*
 * returns 1 if node1 is adjacent to node2, else 0.
 */
int areAdjacent(int node1ID, int node2ID, Node nodes[]) {
  Node *node1 = &nodes[node1ID];
  for (int i = 0; i < node1->nNeighbors; ++i) {
    if (node1->neighbors[i] == node2ID)
      return 1;
  }
  return 0;
}

void addRandomEdge(const int nNodes, Node nodes[]) {
  // iterate in a cycle through all the nodes until an edge is added,
  // starting at random nodes.
  int startNode1 = rand() % nNodes;
  int startNode2 = rand() % nNodes;
  for (int node1 = startNode1;; node1 = (node1 + 1) % nNodes) {
    int tries = 0;
    for (int node2 = startNode2;; node2 = (node2 + 1) % nNodes) {
      tries++;
      // if node 1 is adjacent to all other nodes
      if (tries == nNodes)
        break;
      // don't create an edge to and from the same node
      if (node1 == node2)
        continue;
      if (areAdjacent(node1, node2, nodes))
        continue;

      Node *n1 = &nodes[node1];
      Node *n2 = &nodes[node2];
      //add node 2 to the neighbors of node 1
      n1->neighbors[n1->nNeighbors++] = n2->id;
      //add node 1 to the neighbors of node 2
      n2->neighbors[n2->nNeighbors++] = n1->id;
      return;
    }
  }
}

void assign(Node nodes[], const int conflictingIdSet[], int conflictingSetSize) {
#pragma omp parallel for
  for (int i = 0; i < conflictingSetSize; ++i) {
    Node *conflictingNode = &nodes[conflictingIdSet[i]];
    int smallestAdjacentColor = 1;
    int foundSmallest;
    // keep incrementing smallestAdjacentColor until we find one that isn't used by a neighbor (at the time of reading)
    do {
      foundSmallest = 1;
      for (int j = 0; j < conflictingNode->nNeighbors; ++j) {
        int neighborID = conflictingNode->neighbors[j];
#pragma omp atomic read
        int neighborColor = nodes[neighborID].color;
        if (neighborColor == smallestAdjacentColor) {
          foundSmallest = 0;
          smallestAdjacentColor++;
          break;
        }
      }
    } while (!foundSmallest);
#pragma omp atomic write
    conflictingNode->color = smallestAdjacentColor;
  }
}

int *detectConflicts(int *newConflictingSetSize, Node nodes[], const int conflictingIdSet[], int conflictingSetSize) {
  int *newConflictingIdSet = malloc(sizeof(*newConflictingIdSet) * conflictingSetSize);
  int theNewConflictingSetSize = 0;

  for (int i = 0; i < conflictingSetSize; ++i) {
    Node *conflictingNode = &nodes[conflictingIdSet[i]];
    for (int j = 0; j < conflictingNode->nNeighbors; ++j) {
      Node *neighbor = &nodes[conflictingNode->neighbors[j]];
      if (conflictingNode->color == neighbor->color && conflictingNode->id > neighbor->id) {
        int insertIdx = theNewConflictingSetSize++;
        newConflictingIdSet[insertIdx] = conflictingNode->id;
      }
    }
  }

  *newConflictingSetSize = theNewConflictingSetSize;
  return newConflictingIdSet;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    perror("Program must be ran like ./a.out Num_Nodes Num_Edges Num_Threads");
    exit(1);
  }

  const int nNodes = atoi(argv[1]);
  const int nEdges = atoi(argv[2]);
  const int nThreads = atoi(argv[3]);

  // create a node lookup table.
  // Nodes with id=x will be at index x
  Node *nodes = malloc(sizeof(*nodes) * nNodes);

  for (int i = 0; i < nNodes; ++i) {
    // keeping things simple and allocating memory nNodes-1 number of neighbors
    nodes[i] = (Node) {i, 0, 0, malloc(sizeof(int) * (nNodes - 1))};
  }

  // add random edges
  for (int i = 0; i < nEdges; ++i) {
    addRandomEdge(nNodes, nodes);
  }

  int conflictingSetSize = nNodes;
  int *conflictingIdSet = malloc(sizeof(*conflictingIdSet) * conflictingSetSize);
  for (int i = 0; i < nNodes; ++i) {
    conflictingIdSet[i] = i;
  }

  puts("Starting algo");
  struct timeb start, end;
  ftime(&start);
  // algorithm start
  while (conflictingSetSize != 0) {
    assign(nodes, conflictingIdSet, conflictingSetSize);
    int *newConflictingIdSet = detectConflicts(&conflictingSetSize, nodes, conflictingIdSet, conflictingSetSize);
    free(conflictingIdSet);
    conflictingIdSet = newConflictingIdSet;
  }
  ftime(&end);

  long long diff = (1000 * (end.time - start.time) + (end.millitm - start.millitm));

  printf("%lld", diff);

  //free node neighbors
  for (int i = 0; i < nNodes; ++i) {
    free(nodes[i].neighbors);
  }
  free(nodes);

  return 0;
}