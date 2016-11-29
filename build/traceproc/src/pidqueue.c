#include "pidqueue.h"

#ifndef TRACEPROC_PIDQUEUE_C_
#define TRACEPROC_PIDQUEUE_C_


/*
 * Queue functions.
 */

/*
 * Append a pid to a queue.
 */
void appendPidQueue(PidQueue *pidQueue, pid_t pid) {

	/* Declarations */
	PidQueueNode *pidNode;

	/* Init */
	pidNode = malloc(sizeof(PidQueueNode));
	pidNode->pid = pid;
	pidNode->next = NULL;

	/* Append */
	if (pidQueue->head == NULL) {
		pidQueue->head = pidQueue->tail = pidNode;

	} else {
		pidQueue->tail->next = pidNode;
		pidQueue->tail = pidNode;
	}

	return;
}

/*
 * Get a pid from this queue. If the queue is empty, return -1.
 */
pid_t takePidQueue(PidQueue *pidQueue) {

	/* Declarations */
	pid_t pid;
	PidQueueNode *pidNode;

	/* Check for an empty queue. */
	if (pidQueue->head == NULL)
		return -1;

	/* Get pid */
	pidNode = pidQueue->head;
	pid = pidNode->pid;

	/* Clean up */
	pidQueue->head = pidNode->next;
	free(pidNode);

	if (pidQueue->head == NULL)
		pidQueue->tail = NULL;

	/* Return pid */
	return pid;
}

/*
 * Create a new PidQueue and return a reference to it. This structure should be
 * freed with freePidQueue().
 */
PidQueue *newPidQueue() {
	PidQueue *pidQueue = malloc(sizeof(PidQueue));

	pidQueue->head = pidQueue->tail = NULL;

	return pidQueue;
}

/*
 * Clear queue.
 */
void clearPidQueue(PidQueue *pidQueue) {

	PidQueueNode *nextNode;
	PidQueueNode *thisNode;

	nextNode = pidQueue->head;

	while (nextNode != NULL) {
		thisNode = nextNode;
		nextNode = nextNode->next;

		free(thisNode);
	}

	return;
}

/*
 * Free all resources in a queue.
 */
void freePidQueue(PidQueue *pidQueue) {
	clearPidQueue(pidQueue);
	free(pidQueue);
}

#endif /* TRACEPROC_PIDQUEUE_C_ */
