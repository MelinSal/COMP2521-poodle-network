// COMP2521 Assignment2
//
// proobing a network of computers in a efficient
//
// Authors:
// Melina Salardini (z5393518@unsw.edu.au)
//
// Written: 2/8/2024

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "poodle.h"

// structs

// define a adjacency list
struct AdjListNode {
  int dest;
  int transmissionTime;
  struct AdjListNode *next;
};

struct AdjList {
  struct AdjListNode *head;
};

struct Graph {
  int numComputers;
  struct AdjList *array;
};

// function prototypes
struct probePathResult probePath(struct computer computers[], int numComputers,
                                 struct connection connections[],
                                 int numConnections, int path[],
                                 int pathLength);
static struct AdjListNode *newAdjListNode(int dest, int transmissionTime);
static struct Graph *createGraph(int numComputers);
static void addEdge(struct Graph *graph, int src, int dest,
                    int transmissionTime);
bool *initializeVisited(int numComputers);
static void addPoodleTimeIfNeeded(struct probePathResult *res,
                                  struct computer computers[], bool *visited,
                                  int computer);
static bool hasConnection(struct Graph *graph, int curr, int next);
void addTransmissionTime(struct Graph *graph, int curr, int next,
                         struct probePathResult *res);
static void freeGraph(struct Graph *graph);
struct chooseSourceResult chooseSource(struct computer computers[],
                                       int numComputers,
                                       struct connection connections[],
                                       int numConnections);
static int compare(const void *a, const void *b);
static int BFS(struct Graph *graph, int start, bool *visited,
               int *reachableComputers, struct computer computers[]);
struct poodleResult poodle(struct computer computers[], int numComputers,
                           struct connection connections[], int numConnections,
                           int sourceComputer);
static struct Graph *mst(struct computer computers[], int numComputers,
                         struct connection connections[], int numConnections,
                         int sourceComputer, struct poodleResult *res);
static void initializeArrays(int numComputers, int *isUsedV, int *pathTimes,
                             int sourceComputer, int sourcePoodleTime);
static void initializeFirstStep(struct poodleResult *res, int sourceComputer,
                                int sourcePoodleTime);
static int findCheapestEdge(struct connection connections[], int numConnections,
                            int isUsedV[], struct computer computers[],
                            int pathTimes[], struct poodleResult *res);
static void updateStep(struct poodleResult *res, int computer, int time);
static struct computerList *addRecipient(struct computerList *head,
                                         int computer);
struct poodleResult advancedPoodle(struct computer computers[],
                                   int numComputers,
                                   struct connection connections[],
                                   int numConnections, int sourceComputer);

////////////////////////////////////////////////////////////////////////
// Task 1
// function to return the result of a path we probed
struct probePathResult probePath(struct computer computers[], int numComputers,
                                 struct connection connections[],
                                 int numConnections, int path[],
                                 int pathLength) {
  struct probePathResult res = {SUCCESS, 0};

  // create a graph
  struct Graph *graph = createGraph(numComputers);
  for (int i = 0; i < numConnections; i++) {
    addEdge(graph, connections[i].computerA, connections[i].computerB,
            connections[i].transmissionTime);
  }

  // track visited computers
  bool *visited = initializeVisited(numComputers);

  // iterate through the path
  for (int i = 0; i < pathLength - 1; i++) {
    int curr = path[i];
    int next = path[i + 1];
    if (curr == next) {
      continue;
    }
    addPoodleTimeIfNeeded(&res, computers, visited, curr);

    // check for permission
    if (computers[curr].securityLevel < computers[next].securityLevel - 1) {
      res.status = NO_PERMISSION;
      break;
    }

    if (!hasConnection(graph, curr, next)) {
      res.status = NO_CONNECTION;
      break;
    }

    // Add transmission time
    addTransmissionTime(graph, curr, next, &res);
  }

  // add the last computer if the status is SUCCESS and hasnt been visited
  if (res.status == SUCCESS && !visited[path[pathLength - 1]]) {
    res.elapsedTime += computers[path[pathLength - 1]].poodleTime;
  }

  free(visited);
  freeGraph(graph);
  return res;
}

// a function to create a new adjacency list node
static struct AdjListNode *newAdjListNode(int dest, int transmissionTime) {
  struct AdjListNode *newNode = malloc(sizeof(struct AdjListNode));
  newNode->dest = dest;
  newNode->transmissionTime = transmissionTime;
  newNode->next = NULL;
  return newNode;
}

// a function to create a graph of computers
static struct Graph *createGraph(int numComputers) {
  struct Graph *graph = malloc(sizeof(struct Graph));
  graph->numComputers = numComputers;
  graph->array = malloc(numComputers * sizeof(struct AdjList));
  // initilize the AdjList's head
  for (int i = 0; i < numComputers; i++) {
    graph->array[i].head = NULL;
  }
  return graph;
}

// a function to add edge to a graph
static void addEdge(struct Graph *graph, int src, int dest,
                    int transmissionTime) {
  // add an edge from the source to the destination
  struct AdjListNode *newNode = newAdjListNode(dest, transmissionTime);
  newNode->next = graph->array[src].head;
  graph->array[src].head = newNode;

  // add an edge from destination to source
  newNode = newAdjListNode(src, transmissionTime);
  newNode->next = graph->array[dest].head;
  graph->array[dest].head = newNode;
}

// Function to allocate and initialize the visited array
bool *initializeVisited(int numComputers) {
  bool *visited = calloc(numComputers, sizeof(bool));
  if (visited == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
  }
  return visited;
}

// a function to add the poodle time if needed
static void addPoodleTimeIfNeeded(struct probePathResult *res,
                                  struct computer computers[], bool *visited,
                                  int computer) {
  if (!visited[computer]) {
    res->elapsedTime += computers[computer].poodleTime;
    visited[computer] = true;
  }
}

// a function to check if there is a connection between two nodes
static bool hasConnection(struct Graph *graph, int curr, int next) {
  struct AdjListNode *adjList = graph->array[curr].head;
  while (adjList) {
    if (adjList->dest == next)
      return true;
    adjList = adjList->next;
  }
  return false;
}

// Function to add transmission time
void addTransmissionTime(struct Graph *graph, int curr, int next,
                         struct probePathResult *res) {
  struct AdjListNode *adjList = graph->array[curr].head;
  while (adjList) {
    if (adjList->dest == next) {
      res->elapsedTime += adjList->transmissionTime;
      break;
    }
    adjList = adjList->next;
  }
}

// a function to free the graph
static void freeGraph(struct Graph *graph) {
  for (int i = 0; i < graph->numComputers; i++) {
    struct AdjListNode *current = graph->array[i].head;
    while (current) {
      struct AdjListNode *next = current->next;
      free(current);
      current = next;
    }
  }
  free(graph->array);
  free(graph);
}

////////////////////////////////////////////////////////////////////////
// Task 2
// a function to choose the best source that allows us to
// poodle as many computers as possible
struct chooseSourceResult chooseSource(struct computer computers[],
                                       int numComputers,
                                       struct connection connections[],
                                       int numConnections) {
  struct chooseSourceResult res = {0, 0, NULL};

  // create a graph
  struct Graph *graph = createGraph(numComputers);
  for (int i = 0; i < numConnections; i++) {
    addEdge(graph, connections[i].computerA, connections[i].computerB,
            connections[i].transmissionTime);
  }

  bool *visited = malloc(numComputers * sizeof(bool));
  int *reachableComputers = malloc(numComputers * sizeof(int));

  // for each computer
  for (int i = 0; i < numComputers; i++) {
    // reset the visited array
    for (int j = 0; j < numComputers; j++) {
      visited[j] = false;
    }

    // perform BFS
    int count = BFS(graph, i, visited, reachableComputers, computers);
    // update the source computer if needed
    if (count > res.numComputers) {
      res.sourceComputer = i;
      res.numComputers = count;
      if (res.computers != NULL) {
        free(res.computers);
      }
      // allocate a new array and put the reachable computers inside
      res.computers = malloc(count * sizeof(int));
      for (int j = 0; j < count; j++) {
        res.computers[j] = reachableComputers[j];
      }
      // sort the array of reachable computers
      qsort(res.computers, count, sizeof(int), compare);
    }
  }

  free(visited);
  free(reachableComputers);
  freeGraph(graph);

  return res;
}

// Comparator function for qsort
static int compare(const void *a, const void *b) {
  return (*(int *)a - *(int *)b);
}

// Function to perform BFS and return the number of reachable computers
static int BFS(struct Graph *graph, int start, bool *visited,
               int *reachableComputers, struct computer computers[]) {
  int count = 0;
  int *queue = malloc(graph->numComputers * sizeof(int));
  int front = 0, rear = 0;

  // mark the start node as visited and enqueue the start node
  visited[start] = true;
  queue[rear++] = start;

  // while we havnt explored all the nodes
  while (front < rear) {
    // dequeue the node and mark it as reachable
    int current = queue[front++];
    reachableComputers[count++] = current;

    struct AdjListNode *adjList = graph->array[current].head;
    while (adjList) {
      int adjVertex = adjList->dest;
      if (!visited[adjVertex] && computers[adjVertex].securityLevel - 1 <=
                                     computers[current].securityLevel) {
        visited[adjVertex] = true;
        queue[rear++] = adjVertex;
      }
      adjList = adjList->next;
    }
  }

  free(queue);
  return count;
}

////////////////////////////////////////////////////////////////////////
// Task 3
// a function to poodle as many computers as possible as quick as possible
struct poodleResult poodle(struct computer computers[], int numComputers,
                           struct connection connections[], int numConnections,
                           int sourceComputer) {
  // an array to keep track of which one is poodled
  int *isPoodled = malloc(numComputers * sizeof(int));
  for (int i = 0; i < numComputers; i++) {
    isPoodled[i] = 0;
  }
  isPoodled[sourceComputer] = 1;

  // initilize the res
  struct poodleResult res;
  res.numSteps = 0;
  res.steps = malloc(numComputers * sizeof(struct step));

  // create a mst graph
  struct Graph *mstGraph = mst(computers, numComputers, connections,
                               numConnections, sourceComputer, &res);

  // Update the recipients for each computer by traversing the MST graph
  for (int i = 0; i < res.numSteps; i++) {
    int currentComputer = res.steps[i].computer;
    if (!isPoodled[currentComputer]) {
      continue;
    }

    struct AdjListNode *node = mstGraph->array[currentComputer].head;
    while (node != NULL) {
      if (!isPoodled[node->dest]) {
        res.steps[i].recipients =
            addRecipient(res.steps[i].recipients, node->dest);
        isPoodled[node->dest] = 1;
      }
      node = node->next;
    }
  }

  free(isPoodled);
  freeGraph(mstGraph);
  return res;
}

// a function to create a minimum spanning tree from the given computers
static struct Graph *mst(struct computer computers[], int numComputers,
                         struct connection connections[], int numConnections,
                         int sourceComputer, struct poodleResult *res) {

  struct Graph *mst = createGraph(numComputers);
  int *isUsedV = malloc(numComputers * sizeof(int));
  int *pathTimes = malloc(numComputers * sizeof(int));

  // initilze isUsedV and pathTimes array and the first step
  initializeArrays(numComputers, isUsedV, pathTimes, sourceComputer,
                   computers[sourceComputer].poodleTime);
  initializeFirstStep(res, sourceComputer,
                      computers[sourceComputer].poodleTime);

  while (1) {
    int cheapestEdgeIndex = findCheapestEdge(
        connections, numConnections, isUsedV, computers, pathTimes, res);
    // if no more edge to process break
    if (cheapestEdgeIndex == -1) {
      break;
    }

    // add the edge to mst
    addEdge(mst, connections[cheapestEdgeIndex].computerA,
            connections[cheapestEdgeIndex].computerB,
            connections[cheapestEdgeIndex].transmissionTime);

    // Mark the new computer as used in the isUsedV array
    int computerA;
    int computerB;
    if (!isUsedV[connections[cheapestEdgeIndex].computerA]) {
      computerA = connections[cheapestEdgeIndex].computerB;
      computerB = connections[cheapestEdgeIndex].computerA;
    } else {
      computerA = connections[cheapestEdgeIndex].computerA;
      computerB = connections[cheapestEdgeIndex].computerB;
    }

    int weight = connections[cheapestEdgeIndex].transmissionTime;
    isUsedV[computerB] = 1;
    pathTimes[computerB] =
        pathTimes[computerA] + weight + computers[computerB].poodleTime;
  }

  free(isUsedV);
  free(pathTimes);
  return mst;
}

// a function to initialize the arrays
static void initializeArrays(int numComputers, int *isUsedV, int *pathTimes,
                             int sourceComputer, int sourcePoodleTime) {
  for (int i = 0; i < numComputers; i++) {
    isUsedV[i] = 0;
    pathTimes[i] = INT_MAX;
  }
  isUsedV[sourceComputer] = 1;
  pathTimes[sourceComputer] = sourcePoodleTime;
}

// a function to initilze the first step
static void initializeFirstStep(struct poodleResult *res, int sourceComputer,
                                int sourcePoodleTime) {
  res->steps[res->numSteps].computer = sourceComputer;
  res->steps[res->numSteps].time = sourcePoodleTime;
  res->steps[res->numSteps].recipients = NULL;
  res->numSteps++;
}

// a function to choose the cheapest edge
static int findCheapestEdge(struct connection connections[], int numConnections,
                            int isUsedV[], struct computer computers[],
                            int pathTimes[], struct poodleResult *res) {

  int minTime = INT_MAX;
  int minEdgeIndex = -1;
  int minComputerB;

  for (int i = 0; i < numConnections; i++) {
    int computerA;
    int computerB;
    if (isUsedV[connections[i].computerA] &&
        !isUsedV[connections[i].computerB]) {
      computerA = connections[i].computerA;
      computerB = connections[i].computerB;
    } else if (!isUsedV[connections[i].computerA] &&
               isUsedV[connections[i].computerB]) {
      computerA = connections[i].computerB;
      computerB = connections[i].computerA;
    } else {
      continue;
    }

    // Check if the edge can be used based on security levels
    if (computers[computerA].securityLevel >=
        computers[computerB].securityLevel - 1) {
      int pathTime = pathTimes[computerA] + connections[i].transmissionTime +
                     computers[computerB].poodleTime;
      if (pathTime <= minTime) {
        minTime = pathTime;
        minEdgeIndex = i;
        minComputerB = computerB;
      }
    }
  }

  // If we found a minimum edge, add it to the result steps
  if (minEdgeIndex != -1) {
    updateStep(res, minComputerB, minTime);
  }

  return minEdgeIndex;
}

// a function to update the steps
static void updateStep(struct poodleResult *res, int computer, int time) {
  res->steps[res->numSteps].computer = computer;
  res->steps[res->numSteps].time = time;
  res->steps[res->numSteps].recipients = NULL;
  res->numSteps++;
}

// a function to add a recipient to the linked list
static struct computerList *addRecipient(struct computerList *head,
                                         int computer) {
  struct computerList *newNode = malloc(sizeof(struct computerList));
  newNode->computer = computer;
  newNode->next = NULL;

  // If the list is empty or the new computer should be at the head
  if (head == NULL || head->computer >= computer) {
    newNode->next = head;
    return newNode;
  }

  // Find the correct position to insert the new node
  struct computerList *current = head;
  while (current->next != NULL && current->next->computer < computer) {
    current = current->next;
  }

  // Insert the new node
  newNode->next = current->next;
  current->next = newNode;

  return head;
}

////////////////////////////////////////////////////////////////////////
// Task 4

/**
 * Describe your solution in detail here:
 *
 * TODO
 */
struct poodleResult advancedPoodle(struct computer computers[],
                                   int numComputers,
                                   struct connection connections[],
                                   int numConnections, int sourceComputer) {
  struct poodleResult res = {0, NULL};

  return res;
}
