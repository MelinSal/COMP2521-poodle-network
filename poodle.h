// poodle.h

// !!! DO NOT MODIFY THIS FILE !!!

#ifndef POODLE_H
#define POODLE_H

#include <limits.h>

#define INFINITY INT_MAX

#define MAX_SECURITY_LEVEL 10

struct computer {
	int securityLevel;
	int poodleTime;
};

struct connection {
	int computerA;
	int computerB;
	int transmissionTime;
};

////////////////////////////////////////////////////////////////////////
// Task 1

// An enum is basically used to create a set of #defines,
// except they can be collectively be given a type name
typedef enum status {
	SUCCESS,
	NO_CONNECTION,
	NO_PERMISSION,
} Status;

struct probePathResult {
	Status status;
	int elapsedTime;
};

struct probePathResult probePath(
	struct computer computers[], int numComputers,
	struct connection connections[], int numConnections,
	int path[], int pathLength
);

////////////////////////////////////////////////////////////////////////
// Task 2

struct chooseSourceResult {
	int sourceComputer;
	int numComputers;
	int *computers;
};

struct chooseSourceResult chooseSource(
	struct computer computers[], int numComputers,
	struct connection connections[], int numConnections
);

////////////////////////////////////////////////////////////////////////
// Task 3 and 4

struct poodleResult {
	int numSteps;
	struct step *steps;
};

struct step {
	int computer;
	int time;
	struct computerList *recipients; // only required in task 3
};

struct computerList {
	int computer;
	struct computerList *next;
};

// Task 3
struct poodleResult poodle(
	struct computer computers[], int numComputers,
	struct connection connections[], int numConnections,
	int startingComputer
);

// Task 4
struct poodleResult advancedPoodle(
	struct computer computers[], int numComputers,
	struct connection connections[], int numConnections,
	int startingComputer
);

////////////////////////////////////////////////////////////////////////

#endif

