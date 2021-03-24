#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "sys/time.h"

typedef struct Node {
    // the id of this node (used for total ordering of nodes)
    int id;
    // color of the node
    int color;
    // the number of neighbors
    int nNeighbors;
    // neighbors is an array of Node pointers
    struct Node **neighbors;
} Node;


/*
 * returns 1 if node1 is adjacent to node2, else 0.
 */
int areAdjacent(Node *node1, Node *node2) {
  for (int i = 0; i < node1->nNeighbors; ++i) {
    if (node1->neighbors[i] == node2)
      return 1; //true
  }
  return 0; // false
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

      Node *n1 = &nodes[node1];
      Node *n2 = &nodes[node2];
      if (areAdjacent(n1, n2))
        continue;

      //add node 2 to the neighbors of node 1
      n1->neighbors[n1->nNeighbors++] = n2;
      //add node 1 to the neighbors of node 2
      n2->neighbors[n2->nNeighbors++] = n1;
      return;
    }
  }
}

void assign(Node *conflictingSet[], int conflictingSetSize) {
  #pragma omp parallel for
  for (int i = 0; i < conflictingSetSize; ++i) {
    Node *conflictingNode = conflictingSet[i];
    int smallestAdjacentColor = 1;
    int foundSmallest;
    // keep incrementing smallestAdjacentColor until we find one that isn't used by a neighbor (at the time of reading)
    do {
      foundSmallest = 1; //true
      for (int j = 0; j < conflictingNode->nNeighbors; ++j) {
        Node *neighbor = conflictingNode->neighbors[j];
        int neighborColor;
        #pragma omp atomic read
        neighborColor = neighbor->color;
        if (neighborColor == smallestAdjacentColor) {
          foundSmallest = 0; // false
          smallestAdjacentColor++;
          break;
        }
      }
    } while (!foundSmallest);
    #pragma omp atomic write
    conflictingNode->color = smallestAdjacentColor;
  }
}

Node **detectConflicts(int *newConflictingSetSize, Node *conflictingSet[], int conflictingSetSize) {
  Node **newConflictingSet = malloc(sizeof(Node*) * conflictingSetSize);
  int theNewConflictingSetSize = 0;

  #pragma omp parallel for
  for (int i = 0; i < conflictingSetSize; ++i) {
    Node *conflictingNode = conflictingSet[i];
    for (int j = 0; j < conflictingNode->nNeighbors; ++j) {
      Node *neighbor = conflictingNode->neighbors[j];
      if (conflictingNode->color == neighbor->color && conflictingNode->id > neighbor->id) {
        int insertIdx;
        #pragma omp atomic capture
        insertIdx = theNewConflictingSetSize++;
        newConflictingSet[insertIdx] = conflictingNode;
      }
    }
  }

  *newConflictingSetSize = theNewConflictingSetSize;
  return newConflictingSet;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    perror("Program must be ran like ./a.out Num_Nodes Num_Edges Num_Threads");
    exit(1);
  }

  const int nNodes = atoi(argv[1]);
  const int nEdges = atoi(argv[2]);
  const int nThreads = atoi(argv[3]);

  //checks
  if (nNodes < 4) {
    perror("Program requires at least 4 nodes");
    exit(2);
  }
  if (nEdges > (long)nNodes * (nNodes - 1) / 2) {
    perror("Too many edges");
    exit(3);
  }
  if (nThreads < 1) {
    perror("Program requires at least 1 thread");
    exit(4);
  }

  omp_set_dynamic(0); //don't use dynamic threads
  omp_set_num_threads(nThreads);  //always use nThreads

  // allocate all the nodes
  Node *nodes = malloc(sizeof(Node) * nNodes);
  // initialize node data: id, color, number of neighbors, allocate memory for neighbors
  for (int i = 0; i < nNodes; ++i) {
    // keeping things simple and allocating memory nNodes-1 number of neighbors, this avoid reallocating memory or
    // using a linked list
    nodes[i] = (Node) {i, 0, 0, malloc(sizeof(Node *) * (nNodes - 1))};
  }

  // add random edges
  puts("Building random graph");
  for (int i = 0; i < nEdges; ++i) {
    addRandomEdge(nNodes, nodes);
  }

  int conflictingSetSize = nNodes;
  Node **conflictingSet = malloc(sizeof(Node*) * conflictingSetSize);
  for (int i = 0; i < nNodes; ++i) {
    conflictingSet[i] = &nodes[i];
  }

  puts("Coloring graph");
  struct timeval start, end;
  gettimeofday(&start, NULL);
  // coloring start
  while (conflictingSetSize != 0) {
    assign(conflictingSet, conflictingSetSize);
    Node **newConflictingSet = detectConflicts(&conflictingSetSize, conflictingSet, conflictingSetSize);
    free(conflictingSet);
    conflictingSet = newConflictingSet;
  }
  // coloring end
  gettimeofday(&end, NULL);
  puts("Done coloring");

  long diff = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;

  printf("Time taken: %ld milliseconds\n", diff);

  //print out maximum node degree and maximum color used
  int maxNodeDegree = 0;
  int maxColorUsed = 0;
  for (int i = 0; i < nNodes; ++i) {
    int nodeDegree = nodes[i].nNeighbors;
    int nodeColor = nodes[i].color;
    if (nodeDegree > maxNodeDegree)
      maxNodeDegree = nodeDegree;
    if (nodeColor > maxColorUsed)
      maxColorUsed = nodeColor;
  }
  printf("Max node degree: %d\n", maxNodeDegree);
  printf("Max color used: %d\n", maxColorUsed);

  //de-allocate node neighbors
  for (int i = 0; i < nNodes; ++i) {
    free(nodes[i].neighbors);
  }
  //de-allocate all nodes
  free(nodes);

  return 0;
}