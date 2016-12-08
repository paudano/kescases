#ifndef TRACEPROC_PIDQUEUE_H_
#define TRACEPROC_PIDQUEUE_H_

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

/* An element in a list of pids. */
typedef struct PidQueueNode {
	pid_t pid;
	struct PidQueueNode *next;
} PidQueueNode;

typedef struct PidQueue {
	PidQueueNode *head;
	PidQueueNode *tail;
} PidQueue;

/** Queue functions. */
void appendPidQueue(PidQueue *pidQueue, pid_t pid);
pid_t takePidQueue(PidQueue *pidQueue);
PidQueue *newPidQueue();
void clearPidQueue(PidQueue *pidQueue);
void freePidQueue(PidQueue *pidQueue);

#endif /* TRACEPROC_PIDQUEUE_H_ */
