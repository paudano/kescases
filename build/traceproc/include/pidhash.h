#ifndef TRACEPROC_PIDHASH_H_
#define TRACEPROC_PIDHASH_H_

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TRACEPROC_PIDHASH_REPROBE_MAX 8

/*
 * pid_node_t: Information about a process in a list of pids with the same
 * parent pid.
 */
typedef struct PidNode {  /* http://man7.org/linux/man-pages/man5/proc.5.html */
	pid_t pid;             /* Process ID */
	pid_t ppid;            /* Parent process ID */
	char state;            /* Process state (R, S, D, Z, T, etc) */
	unsigned long vsize;   /* Virtual size of process (bytes) */
	long rss;              /* Resident size (pages) */
	struct PidNode *next;  /* Next node in a linked-list of PidNodes */
} PidNode;

/*
 * A list of child pids of a given parent.
 */
typedef struct ParentPidNode {
	pid_t ppid;      /* Parent process ID */
	PidNode *child;  /* Head of a linked-list of child nodes associated with this process */
} ParentPidNode;

/*
 * Structures needed for hashing.
 */
typedef struct PidHash {
	ParentPidNode * const hashArray;  /* Array of parent nodes */
	char * const nameBuffer;          /* Buffer for storing file names in /proc */
	char * const statBuffer;          /* Buffer for reading /proc/[pid]/stat (must be large enough to hold the whole file) */
	const int hashArraySize;          /* Size of hashArray */
	const int nameBufferSize;         /* Size of nameBuffer */
	const int statBufferSize;         /* Size of statBuffer */
} PidHash;

/* Hash routines */
#define traceproc_pidhash_hash(PID) (abs((int) ((long) (PID) * 2654435761)))
#define traceproc_pidhash_reprobe(PID) (((PID) << 7 ^ (PID) << 2 ^ (PID) << 1) >> 1)

/*
 * Function prototypes.
 */

PidHash *newPidHash(int hashArraySize, int nameBufferSize, int statBufferSize);
void freePidHash(PidHash *pidHashPtr);
void clearPidHash(PidHash *pidHashPtr);

pid_t addPidHash(PidHash *pidHashPtr, pid_t pid);
PidNode *getPidNode(PidHash *pidHashPtr, pid_t pid);
int getPidHashIndex(PidHash *pidHashPtr, pid_t pid);

#endif /* TRACEPROC_PIDHASH_H_ */
