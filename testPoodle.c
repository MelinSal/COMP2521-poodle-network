// Main program for testing

// !!! DO NOT MODIFY THIS FILE !!!

#include <assert.h>
#include <err.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "poodle.h"

#define MAX_PRINTABLE_PROBE_LEN 50

struct test {
	int numComputers;
	struct computer *computers;
	int numConnections;
	struct connection *connections;
	int pathLength;
	int *path;
	int source;
};

static void testProbePath(FILE *fp);
static void testChooseSource(FILE *fp);
static void testPoodle(FILE *fp);
static void testAdvancedPoodle(FILE *fp);

static void scanNetwork(FILE *fp, struct test *t);
static void readNetworkFile(char *filename, struct test *t);

static void scanProbePath(FILE *fp, struct test *t);
static void scanSource(FILE *fp, struct test *t);

int main(int argc, char *argv[]) {
	if (argc != 2) {
		errx(EXIT_FAILURE, "usage: %s <test file>", argv[0]);
	}

	FILE *fp = fopen(argv[1], "r");
	if (fp == NULL) {
		errx(EXIT_FAILURE, "error: failed to open '%s' for reading\n",
		     argv[1]);
	}

	int task = 0;
	if (fscanf(fp, "%d", &task) != 1) {
		errx(EXIT_FAILURE, "error: failed to read task number\n");
	}

	switch (task) {
		case 1:
			testProbePath(fp);
			break;
		case 2:
			testChooseSource(fp);
			break;
		case 3:
			testPoodle(fp);
			break;
		case 4:
			testAdvancedPoodle(fp);
			break;
		default:
			errx(EXIT_FAILURE, "error: invalid task number '%d'", task);
	}

	fclose(fp);
}

static void testProbePath(FILE *fp) {
	struct test t = {};
	scanNetwork(fp, &t);
	scanProbePath(fp, &t);

	printf("Probe path: ");
	if (t.pathLength <= MAX_PRINTABLE_PROBE_LEN) {
		for (int i = 0; i < t.pathLength; i++) {
			if (i > 0) {
				printf(" -> ");
			}
			printf("%d", t.path[i]);
		}
	} else {
		printf("[too long to print]");
	}
	printf("\n\n");

	struct probePathResult res = probePath(
		t.computers, t.numComputers,
		t.connections, t.numConnections,
		t.path, t.pathLength
	);

	printf("Result: ");
	if (res.status == SUCCESS) {
		printf("Success\n");
	} else {
		printf("Error\n");
		printf("Reason: ");
		if (res.status == NO_CONNECTION) {
			printf("No connection\n");
		} else if (res.status == NO_PERMISSION) {
			printf("No permission\n");
		} else {
			printf("Unknown reason\n");
		}
	}
	printf("Elapsed time: %d seconds\n", res.elapsedTime);

	free(t.path);
	free(t.connections);
	free(t.computers);
}

static void testChooseSource(FILE *fp) {
	struct test t = {};
	scanNetwork(fp, &t);

	struct chooseSourceResult res = chooseSource(
		t.computers, t.numComputers,
		t.connections, t.numConnections
	);

	printf("Source computer: computer %d\n", res.sourceComputer);
	printf("Number of reachable computers: %d\n", res.numComputers);
	printf("Reachable computers:\n");
	for (int i = 0; i < res.numComputers; i++) {
		printf("- computer %d\n", res.computers[i]);
	}

	free(res.computers);

	free(t.connections);
	free(t.computers);
}

static void testPoodle(FILE *fp) {
	struct test t = {};
	scanNetwork(fp, &t);
	scanSource(fp, &t);

	struct poodleResult res = poodle(
		t.computers, t.numComputers,
		t.connections, t.numConnections,
		t.source
	);

	printf("Plan:\n");
	for (int i = 0; i < res.numSteps; i++) {
		printf("- computer %d poodled at %d seconds\n",
		       res.steps[i].computer, res.steps[i].time);
		
		if (res.steps[i].recipients != NULL) {
			printf("  - pug sent to: ");
			struct computerList *curr = res.steps[i].recipients;
			for (; curr != NULL; curr = curr->next) {
				printf("%d", curr->computer);
				if (curr->next != NULL) {
					printf(", ");
				}
			}
			printf("\n");
		}
	}

	for (int i = 0; i < res.numSteps; i++) {
		struct computerList *curr = res.steps[i].recipients;
		while (curr != NULL) {
			struct computerList *temp = curr;
			curr = curr->next;
			free(temp);
		}
	}
	free(res.steps);

	free(t.connections);
	free(t.computers);
}

static void testAdvancedPoodle(FILE *fp) {
	struct test t = {};
	scanNetwork(fp, &t);
	scanSource(fp, &t);

	struct poodleResult res = advancedPoodle(
		t.computers, t.numComputers,
		t.connections, t.numConnections,
		t.source
	);

	printf("Plan:\n");
	for (int i = 0; i < res.numSteps; i++) {
		printf("- computer %d poodled at %d seconds\n",
		       res.steps[i].computer, res.steps[i].time);
	}
	free(res.steps);

	free(t.connections);
	free(t.computers);
}

////////////////////////////////////////////////////////////////////////

static void scanNetwork(FILE *fp, struct test *t) {
	char filename[100] = {0};
    if (fscanf(fp, "%s", filename) != 1) {
        errx(EXIT_FAILURE, "error: failed to read network filename");
    }
    char path[100] = {0};
    sprintf(path, "data/%s", filename);

	readNetworkFile(path, t);
}

static void readNetworkFile(char *filename, struct test *t) {
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		errx(EXIT_FAILURE, "error: failed to open '%s' for reading",
		     filename);
	}

	if (fscanf(fp, "%d", &(t->numComputers)) != 1) {
		errx(EXIT_FAILURE, "error: failed to read number of computers");
	}

	if (fscanf(fp, "%d", &(t->numConnections)) != 1) {
		errx(EXIT_FAILURE, "error: failed to read number of connections");
	}

	t->computers = malloc(t->numComputers * sizeof(struct computer));
	if (t->computers == NULL) {
		errx(EXIT_FAILURE, "error: out of memory");
	}

	for (int i = 0; i < t->numComputers; i++) {
		int securityLevel;
		int poodleTime;
		
		if (fscanf(fp, "%d %d", &securityLevel, &poodleTime) != 2) {
			errx(EXIT_FAILURE, "error: failed to read computer details");
		} else if (securityLevel < 1 || securityLevel > MAX_SECURITY_LEVEL) {
			errx(EXIT_FAILURE, "error: invalid security level '%d'",
			     securityLevel);
		} else if (poodleTime <= 0) {
			errx(EXIT_FAILURE, "error: invalid poodle time '%d'", poodleTime);
		}
		
		t->computers[i] = (struct computer){securityLevel, poodleTime};
	}

	t->connections = malloc(t->numConnections * sizeof(struct connection));
	if (t->connections == NULL) {
		errx(EXIT_FAILURE, "error: out of memory");
	}

	for (int i = 0; i < t->numConnections; i++) {
		int compA;
		int compB;
		int transmissionTime;

		if (fscanf(fp, "%d %d %d", &compA, &compB, &transmissionTime) != 3) {
			errx(EXIT_FAILURE, "error: failed to read connection details");
		} else if (compA < 0 || compA >= t->numComputers) {
			errx(EXIT_FAILURE, "error: invalid computer number '%d'", compA);
		} else if (compB < 0 || compB >= t->numComputers) {
			errx(EXIT_FAILURE, "error: invalid computer number '%d'", compB);
		} else if (transmissionTime <= 0) {
			errx(EXIT_FAILURE, "error: invalid transmission time '%d'",
			     transmissionTime);
		}

		t->connections[i] = (struct connection){compA, compB, transmissionTime};
	}
	
	fclose(fp);
}

////////////////////////////////////////////////////////////////////////

static void scanProbePath(FILE *fp, struct test *t) {
	if (fscanf(fp, "%d", &(t->pathLength)) != 1) {
		errx(EXIT_FAILURE, "error: failed to read path length");
	}

	t->path = malloc(t->pathLength * sizeof(int));
	if (t->path == NULL) {
		errx(EXIT_FAILURE, "error: out of memory");
	}

	for (int i = 0; i < t->pathLength; i++) {
		if (fscanf(fp, "%d", &(t->path[i])) != 1) {
			errx(EXIT_FAILURE, "error: failed to read probe path");
		}
	}
}

static void scanSource(FILE *fp, struct test *t) {
	if (fscanf(fp, "%d", &(t->source)) != 1) {
		errx(EXIT_FAILURE, "error: failed to read source computer");
	}
}

