#include <stdio.h>
#include <errno.h>

#include <pidhash.h>


/*
 * Create a PidHash structure.
 */
PidHash *newPidHash(int hashArraySize, int nameBufferSize, int statBufferSize) {

	/* Buffers */
	ParentPidNode *hashArray;
	char * nameBuffer;
	char * statBuffer;
	PidHash *pidHashPtr;

	/* Check arguments */
	if (hashArraySize < 0 || nameBufferSize < 0 || statBufferSize < 0) {
		errno = EINVAL;
		return NULL;
	}

	/* nameBuffer init */
	nameBuffer = calloc(nameBufferSize, sizeof(char));

	if (nameBuffer == NULL)
		return NULL;

	/* statBuffer init */
	statBuffer = calloc(statBufferSize, sizeof(char));

	if (statBuffer == NULL) {
		free(nameBuffer);

		return NULL;
	}

	/* hashArray allocation */
	hashArray = calloc(hashArraySize, sizeof(ParentPidNode));

	if (hashArray == NULL) {
		free(nameBuffer);
		free(statBuffer);

		return NULL;
	}

	/* Initialize array. */
	for (int index = 0; index < hashArraySize; ++index) {
		hashArray[index].ppid = -1;
		hashArray[index].child = NULL;
	}

	/* Memory for PidHash structure */
	pidHashPtr = malloc(sizeof(PidHash));

	if (pidHashPtr == NULL) {
		free(nameBuffer);
		free(statBuffer);
		free(hashArray);

		return NULL;
	}

	/* Initialize structure */
	PidHash pidHash = {
			.hashArray = hashArray, .hashArraySize = hashArraySize,
			.nameBuffer = nameBuffer, .nameBufferSize = nameBufferSize,
			.statBuffer = statBuffer, .statBufferSize = statBufferSize
	};

	/* Copy to heap memory. */
	memcpy(pidHashPtr, &pidHash, sizeof(PidHash));

	/* Return object on the heap. */
	return pidHashPtr;
}

/*
 * Free a PidHash structure and all its resources.
 */
void freePidHash(PidHash *pidHashPtr) {

	if (pidHashPtr == NULL)
		return;

	clearPidHash(pidHashPtr);

	free(pidHashPtr->nameBuffer);
	free(pidHashPtr->statBuffer);
	free(pidHashPtr->hashArray);
	free(pidHashPtr);

	return;
}

/*
 * Clear a PidHash of all child elements.
 */
void clearPidHash(PidHash *pidHashPtr) {

	PidNode *pidNodePtr;
	PidNode *pidNodePtrToFree;

	if (pidHashPtr == NULL)
		return;

	/* Traverse hash. */
	for (int index = 0; index < pidHashPtr->hashArraySize; ++index) {

		/* Traverse child chain. */
		pidNodePtr = pidHashPtr->hashArray[index].child;

		while (pidNodePtr != NULL) {
			pidNodePtrToFree = pidNodePtr;
			pidNodePtr = pidNodePtr->next;

			free(pidNodePtrToFree);
		}

		pidHashPtr->hashArray[index].ppid = -1;
		pidHashPtr->hashArray[index].child = NULL;
	}
}

/*
 * Add a pid to the hash and return the parent pid. If pid is 0, -1 is
 * returned.
 *
 * Args:
 *   * pid: Process ID.
 *   * pidHashPtr: A pointer to a structure of PidHash (buffers and hash).
 *
 * Return: Parent pid of pid, -1 if pid is 0, or -2 if an error occurs (errno set).
 *
 * -2 is returned on error and errno is set:
 *   * EINVAL: pid is less than 0 or pidHashPtr is NULL.
 *   * EOVERFLOW: PID hash is full.
 *   * Other: As set by other calls (such as fopen).
 */
pid_t addPidHash(PidHash *pidHashPtr, pid_t pid) {

	PidNode *pidNodePtr;  /* Pointer to a PidNode structure to be added to the array. */

	int hashIndex;     /* Index of the hash array where the element is found. */
	int reprobeCount;  /* Number of reprobe attempts. */
	int pidKey;        /* Key into the hash array. This is initially the pid, but is altered when reprobing. */

	/* Check arguments. */
	if (pid < 0 || pidHashPtr == NULL) {
		errno = EINVAL;
		return -2;
	}

	/* Get pid information. */
	pidNodePtr = getPidNode(pidHashPtr, pid);

	if (pidNodePtr == NULL)
		return -2;  // errno set

	/* Find parent element or an empty space in the hash array. */
	pidKey = pidNodePtr->ppid;
	hashIndex = traceproc_pidhash_hash(pidKey) % pidHashPtr->hashArraySize;
	reprobeCount = 0;

	while (pidHashPtr->hashArray[hashIndex].ppid != pidNodePtr->ppid &&
		   pidHashPtr->hashArray[hashIndex].ppid >= 0 &&
		   reprobeCount < TRACEPROC_PIDHASH_REPROBE_MAX) {

		pidKey = traceproc_pidhash_reprobe(pidKey);
		hashIndex = traceproc_pidhash_hash(pidKey) % pidHashPtr->hashArraySize;
		++reprobeCount;
	}

	/* Parent pid not found, add. */
	if (pidHashPtr->hashArray[hashIndex].ppid < 0) {
		pidHashPtr->hashArray[hashIndex].ppid = pidNodePtr->ppid;
		pidHashPtr->hashArray[hashIndex].child = pidNodePtr;

		if (pid == 0)
			return -1;

		return pidNodePtr->ppid;
	}

	/* Ran out of reprobe attempts. */
	if (pidHashPtr->hashArray[hashIndex].ppid != pidNodePtr->ppid) {
		errno = EOVERFLOW;
		return -2;
	}

	/* Parent pid node found, add pid. */
	pidNodePtr->next = pidHashPtr->hashArray[hashIndex].child;
	pidHashPtr->hashArray[hashIndex].child = pidNodePtr;

	if (pid == 0)
		return -1;

	return pidNodePtr->ppid;
}

/*
 * Get index of pidHashPtr->hashArray where ppid is located.
 *
 * Return: Index of pidHashPtr->hashArray, -1 if ppid was not found, and -2 on error.
 */
int getPidHashIndex(PidHash *pidHashPtr, pid_t ppid) {

	int hashIndex;
	int pidKey;
	int reprobeCount;

	/* Check arguments. */
	if (ppid < 0 || pidHashPtr == NULL) {
		errno = EINVAL;
		return -2;
	}

	/* Find parent element or an empty space in the hash array. */
	pidKey = ppid;
	hashIndex = traceproc_pidhash_hash(pidKey) % pidHashPtr->hashArraySize;
	reprobeCount = 0;

	while (pidHashPtr->hashArray[hashIndex].ppid != ppid &&
		   pidHashPtr->hashArray[hashIndex].ppid >= 0 &&
		   reprobeCount < TRACEPROC_PIDHASH_REPROBE_MAX) {

		pidKey = traceproc_pidhash_reprobe(pidKey);
		hashIndex = traceproc_pidhash_hash(pidKey) % pidHashPtr->hashArraySize;
		++reprobeCount;
	}

	/* ppid not found. */
	if (pidHashPtr->hashArray[hashIndex].ppid != ppid)
		return -1;

	/* Found. */
	return hashIndex;
}

/*
 * Get a node of PID information. The structure returned must be released with a call to
 * free().
 *
 * Args:
 *   * pid: Process ID.
 *   * pidHashPtr: A pointer to a structure of PidHash (buffers and hash).
 *
 * NULL is returned on error and errno is set:
 *   * EINVAL: pid is less than 0 or pidHashPtr is NULL.
 *   * EOVERFLOW: PID hash is full.
 *   * EIO: I/O error reading /proc/[pid]/stat
 *   * Other: As set by other calls (such as fopen).
 */
PidNode *getPidNode(PidHash *pidHashPtr, pid_t pid) {

	int ppid;             /* Parent PID. */
	char state;           /* Process state. */
	unsigned long vsize;  /* Virtual size (bytes). */
	long rss;             /* RSS (bytes). */

	FILE *statFile;  /* Status file (/proc/PID/stat) for the child process. */

	int nRead;  /* Number of bytes read into statBuffer. */
	int nScan;  /* Number of bytes interpreted by sscanf. */

	int field;  /* Field number currently being scanned. */

	char *statPtr;

	PidNode *pidNode;  /* Node to return. */

	/* Init */
	field = 1;

	/* Open stat. */
	snprintf(pidHashPtr->nameBuffer, pidHashPtr->nameBufferSize, "/proc/%d/stat", pid);

	statFile = fopen(pidHashPtr->nameBuffer, "r");

	if (statFile == NULL)
		return NULL;

	/* Read stat. */
	nRead = fread(pidHashPtr->statBuffer, sizeof(char), pidHashPtr->statBufferSize, statFile);

	if (nRead < 0)
		return NULL;

	statPtr = pidHashPtr->statBuffer;

	while (field < 25) {

		/* command (may have whitepsace in the name). */
		if (field == 2) {
			while (*statPtr != ')')
				++statPtr;

			while (*statPtr != ' ')
				++statPtr;

			++statPtr;
			++field;

			continue;
		}

		/* state */
		if (field == 3) {
			state = *statPtr;

			statPtr += 2;
			++field;

			continue;
		}

		/* ppid */
		if (field == 4) {
			nScan = sscanf(statPtr, "%d", &ppid);

			if (nScan == EOF || nScan < 1) {
				errno = EIO;
				return NULL;
			}

			while (*statPtr != ' ')
				++statPtr;

			++statPtr;
			++field;

			continue;
		}

		/* vsize */
		if (field == 23) {
			nScan = sscanf(statPtr, "%lu", &vsize);

			if (nScan == EOF || nScan < 1) {
				errno = EIO;
				return NULL;
			}

			while (*statPtr != ' ')
				++statPtr;

			++statPtr;
			++field;

			continue;
		}

		/* rss */
		if (field == 24) {
			nScan = sscanf(statPtr, "%ld", &rss);

			if (nScan == EOF || nScan < 1) {
				errno = EIO;
				return NULL;
			}

			while (*statPtr != ' ')
				++statPtr;

			++statPtr;
			++field;

			continue;
		}

		/* Other field (skip) */
		while (*statPtr != ' ')
			++statPtr;

		++statPtr;
		++field;
	}

	/* Clean up. */
	fclose(statFile);

	/* Return nothing if the process was not read. */
	if (ppid < 0)
		return NULL;

	/* Create node and return. */
	pidNode = malloc(sizeof(PidNode));

	pidNode->pid = pid;
	pidNode->ppid = ppid;
	pidNode->state = state;
	pidNode->vsize = vsize;
	pidNode->rss = rss;
	pidNode->next = NULL;

	return pidNode;
}
