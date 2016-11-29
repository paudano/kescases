/*
 Track memory at periodic intervals and report memory usage.
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <libgen.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <sys/wait.h>
#include <sched.h>
#include <dirent.h>
#include <ctype.h>

#include "pidhash.h"
#include "pidqueue.h"


/* Constants. */
const int ERR_NONE = 0;
const int ERR_USAGE = 1;
const int ERR_IO = 2;
const int ERR_THREAD = 3;
const int ERR_PARSE = 4;

#define HASH_ARRAY_SIZE 16834
#define NAME_BUF_SIZE   1024
#define STAT_BUF_SIZE   4096

/* Structures */
typedef struct ProcMonitorSummary {
	long runTime;    /* Run time (ms). */
	long maxRss;     /* Max resident size (memory usage) in kB. */
	int exitStatus;  /* Process return code. */
} ProcMonitorSummary;

typedef struct TraceArguments {
	int pid;         /* Process to monitor. */
	FILE *outFile;  /* File time snapshots are written to. */
	long startTime;   /* Start time (ms) since epoch. */
	int waitTime;   /* Time to wait (s) between snapshots. */
	int active;      /* Set to FALSE (zero) to tell the trace process to stop. */
	long maxRss;    /* Return max resident size (memory usage) in kB. */
} TraceArguments;

typedef struct ProcSnapshot {
	long rss;    /* Resident size. */
	long vsize;  /* Virtual size. */
	int nproc;  /* Number of processes spawned (including itself) */
} ProcSnapshot;

typedef struct MonitorThreadStatus {
	int returnValue;  /* Return value of the process (1 for success, 0 for failure) */
	int err;          /* errno set when thread exited. */
	char *msg;        /* Error message if returnValue is 0, otherwise NULL. String is static and should not be freed */
} MonitorThreadStatus;

/* Function prototypes. */
void monitorProcess(FILE *outFile, char **argv, int waitTime, ProcMonitorSummary *monitorSummary, int cpuSet);
int getSnapshot(pid_t pid, ProcSnapshot *procSnapshot, PidQueue *pidQueue, PidHash *pidHash);
void *traceProcess(void *traceArgsVoidPtr);
int isPidDir(const char *charPtr);

void printSummary(ProcMonitorSummary *monitorSummary);

void printHelp(void);

void err(const char *msg, int errCode);
void errnoErr(const char *msg, int errCode);

/* Program name. */
#define PROG_NAME_SIZE 25
char *progName;

/* Page size (set in main()) */
int pageSize;

/* Function: main */
int main(int argc, char **argv) {

	/* Declarations. */
	int option;    /* Option character. */
	int inOpts;    /* Flag to stop option processing when false (0). */
	int waitTime;  /* Number of seconds to wait between checking the child process. */
	int verbose;   /* Flag to output verbose information about the process after it terminates. */

	char *outFileName;  /* Name of the output file. */

	char **endPtr;  /* End pointer for strtol. */

	FILE *outFile;  /* Output file handle */

	ProcMonitorSummary monitorSummary;  /* Summary from process monitoring */

	char progNameLocal[PROG_NAME_SIZE];
	progName = progNameLocal;

	int cpuNum;
	unsigned int cpuSet;

	/* Initialize. */
	inOpts = 1;
	waitTime = 1;
	outFileName = NULL;
	endPtr = NULL;
	verbose = 0;

	cpuSet = 0x00;

	pageSize = getpagesize();

	/* Set program name for error reporting. */
	if (argc > 0)
		strncpy(progName, basename(argv[0]), PROG_NAME_SIZE);
	else
		strncpy(progName, "traceproc", PROG_NAME_SIZE);

	/* Process arguments. */
	while (inOpts == 1 && (option = getopt(argc, argv, "+c:ho:w:v")) != -1) {

		switch (option) {

		case 'h':
			printHelp();
			exit(ERR_NONE);

		case 'w':
			waitTime = strtol(optarg, endPtr, 0);

			if (endPtr != NULL)
				err("Invalid number of seconds to wait between child process checks (-s)", ERR_USAGE);

			if (waitTime < 1)
				err("Wait time must not be less than 1 (-s)", ERR_USAGE);

			break;

		case 'c':
			cpuNum = strtol(optarg, endPtr, 0);

			if (endPtr != NULL)
				err("Invalid CPU number (-c)", ERR_USAGE);

			if (cpuNum < 0)
				err("CPU number must not be less than 0 (-c)", ERR_USAGE);

			if (cpuNum > sizeof(int) * 8)
				err("CPU number is too large (-s)", ERR_USAGE);

			cpuSet = cpuSet | (0x01 << cpuNum);

			break;
		
		case 'v':
			verbose = 1;
			break;

		case 'o':
			outFileName = optarg;
			break;
		}
	}

	if (optind >= argc)
		err("No command to be run (see -h for help)", ERR_USAGE);

	/* Open output file. */
	if (outFileName != NULL) {
		outFile = fopen(outFileName, "w");

		fprintf(outFile, "#n\ttime_ms\trss_kb\tvss_kb\tn_proc\n");

	} else {
		outFile = stdout;
	}

	/* Fork and monitor. */
	monitorProcess(outFile, argv + optind, waitTime, &monitorSummary, cpuSet);

	/* Close file. */
	if (outFile != NULL && outFile != stdout)
		fclose(outFile);

	/* Print summary. */
	if (verbose)
		printSummary(&monitorSummary);

	return monitorSummary.exitStatus;
}

/* Function: monitorProcess */
void monitorProcess(FILE *outFile, char **argv, int waitTime, ProcMonitorSummary *monitorSummary, int cpuSet) {

	/* Declarations. */

	struct timespec currentTime;  /* Timespec struct for getting the current time. */
	long startTime;               /* Clock time: milliseconds at start. */
	long ms;                      /* Elapsed time: milliseconds. */

	int childPid;    /* Monitored PID. */
	int exitStatus;  /* Monitored PID exit status. */

	TraceArguments traceArgs;  /* Arguments for the trace function and write-back data from the trace function. */

	pthread_t traceThread;  /* Thread running the trace */

	cpu_set_t cpuMask;  /* Mask describing the CPUs the process should run on. */

	MonitorThreadStatus *monitorExitStatus;

	/* Fork child process. */
	childPid = fork();

	if (childPid == -1)
		errnoErr("Error creating child process to be monitored", ERR_THREAD);

	/* Child: Exec */
	if (childPid == 0) {

		/* Set CPUS */
		if (cpuSet != 0) {

			CPU_ZERO(&cpuMask);

			for (int index = 0; index < sizeof(unsigned int) * 8; ++index) {
				if ((cpuSet & 0x01 << index) != 0)
					CPU_SET(index, &cpuMask);
			}

			if (sched_setaffinity(0, sizeof(cpuMask), &cpuMask) != 0)
				errnoErr("Error setting CPUs", ERR_THREAD);
		}

		execvp(argv[0], argv);

		errnoErr("Error running process", ERR_THREAD);

		return;
	}

	/* Get start time. */
	clock_gettime(CLOCK_MONOTONIC, &currentTime);

	startTime = currentTime.tv_sec * 1000 + currentTime.tv_nsec / 1000000;

	/* Set traceArgs. */
	traceArgs.pid = childPid;
	traceArgs.outFile = outFile;
	traceArgs.startTime = startTime;
	traceArgs.waitTime = waitTime;
	traceArgs.active = 1;
	traceArgs.maxRss = 0;

	/* Start trace process. */
	if (pthread_create(&traceThread, NULL, &traceProcess, &traceArgs) != 0)
		err("Error creating trace thread", ERR_THREAD);

	/* Wait for process to stop. */
	waitpid(childPid, &exitStatus, 0);

	clock_gettime(CLOCK_MONOTONIC, &currentTime);

	ms = currentTime.tv_sec * 1000 + currentTime.tv_nsec / 1000000;

	/* Wait for trace process to stop. */
	traceArgs.active = 0;

	pthread_join(traceThread, (void **) &monitorExitStatus);

	if (! monitorExitStatus->returnValue) {
		fprintf(outFile, "# Monitor process failed: %s (%s)\n", monitorExitStatus->msg, strerror(monitorExitStatus->err));
	}

	free(monitorExitStatus);

	/* Set monitor info. */
	if (monitorSummary != NULL) {
		monitorSummary->maxRss = traceArgs.maxRss;
		monitorSummary->runTime = ms - startTime;
		monitorSummary->exitStatus = WEXITSTATUS(exitStatus);
	}

	return;
}

/* Function: traceProcess */
void *traceProcess(void *traceArgsVoidPtr) {

	/* Declarations. */
	long count;  /* Number of iterations. */

	long ms;                      /* Elapsed time (milliseconds). */
	struct timespec currentTime;  /* Timespec struct for getting the current time. */

	long sleepMs;                 /* Time to sleep (ms) between querying process usage */
	struct timespec sleepTime;    /* Timespec struct argument nanosleep() */

	pid_t pid;  /* Process ID monitored. */

	TraceArguments *traceArgs;  /* Trace arguments. */

	ProcSnapshot procSnapshot;  /* Snapshot of the process at one point in time. */

	PidQueue *pidQueue;  /* Queue of PIDs (reused by getResidentSize()). */
	PidHash *pidHash;    /* Hash of PIDs organized by parent process (reused by getResidentSize()). */

	MonitorThreadStatus *exitStatus;

	/* Check arguments. */
	if (traceArgsVoidPtr == NULL) {
		pthread_exit(NULL);
	}

	/* Init */
	errno = 0;

	exitStatus = malloc(sizeof(MonitorThreadStatus));
	exitStatus->returnValue = 0;
	exitStatus->err = 0;
	exitStatus->msg = NULL;

	pidQueue = newPidQueue();
	pidHash = newPidHash(HASH_ARRAY_SIZE, NAME_BUF_SIZE, STAT_BUF_SIZE);
	traceArgs = (TraceArguments *) traceArgsVoidPtr;

	if (pidQueue == NULL || pidHash == NULL) {
		exitStatus->returnValue = 0;
		exitStatus->err = errno;

		if (pidQueue == NULL)
			exitStatus->msg = "Failed allocating PID queue";
		else
			exitStatus->msg = "Failed allocating PID hash";

		pthread_exit(exitStatus);
	}

	pid = traceArgs->pid;

	/* Monitor (init). */
	ms = traceArgs->startTime;

	count = 0;

	/* Monitor VSS until the child process dies. */
	while (traceArgs->active) {

		clock_gettime(CLOCK_MONOTONIC, &currentTime);

		/* Get process stats */
		if (! getSnapshot(pid, &procSnapshot, pidQueue, pidHash)) {
			exitStatus->returnValue = 0;
			exitStatus->err = errno;
			exitStatus->msg = "Failed getting process snapshot";

			pthread_exit(exitStatus);
		}

		/* Set max. */
		if (procSnapshot.rss > traceArgs->maxRss)
			traceArgs->maxRss = procSnapshot.rss;

		/* Parse time-stamp. */
		ms = currentTime.tv_sec * 1000 + currentTime.tv_nsec / 1000000 - traceArgs->startTime;

		/* Write. */
		if (traceArgs->outFile != NULL) {
			fprintf(traceArgs->outFile, "%d\t%ld.%ld\t%ld\t%ld\t%d\n",
					++count,
					ms / 1000, ms % 1000,
					procSnapshot.rss * pageSize / 1024,
					procSnapshot.vsize / 1024,
					procSnapshot.nproc
			);

			fflush(traceArgs->outFile);
		}

		/* Sleep to next round. */
		sleepMs = count * traceArgs->waitTime * 1000 - ms;

		if (sleepMs > 10) {

			sleepTime.tv_sec = sleepMs / 1000;
			sleepTime.tv_nsec = sleepMs % 1000 * 1000000;

			/* Pause between cycles. */
			nanosleep(&sleepTime, NULL);
		}
	}

	/* Free resources. */
	freePidQueue(pidQueue);
	freePidHash(pidHash);

	/* Exit. */
	exitStatus->returnValue = 1;

	pthread_exit(exitStatus);
}

/* Function: getResidentSize
 *
 * Get the estimated resident size of the process in memory in kB.
 */
int getSnapshot(pid_t pid, ProcSnapshot *procSnapshot, PidQueue *pidQueue, PidHash *pidHash) {

	/* Declarations. */
	DIR *procDir;               /* Structure for reading /proc */
	struct dirent *procDirEnt;  /* Entry in /proc */

	pid_t pidEntry;  /* Numeric pid in /proc or in the pid queue */

	PidNode *pidNode;  /* Resources of pid */

	int hashIndex;  /* Index of pidHash where a process is found. */

	/* Init */
	clearPidQueue(pidQueue);
	clearPidHash(pidHash);

	procSnapshot->rss = 0;
	procSnapshot->vsize = 0;
	procSnapshot->nproc = 0;

	errno = 0;

	/* Open /proc */
	procDir = opendir("/proc");

	if (procDir == NULL)
		return 0;  // errno set

	/* Read /proc and build pid hash. */
	procDirEnt = readdir(procDir);

	while (procDirEnt != NULL) {
		if (isPidDir(procDirEnt->d_name)) {

			/* Add to hash */
			if (addPidHash(pidHash, atoi(procDirEnt->d_name)) <= -2) {

				// Do not fail if the process is no longer in /proc (died since directory listing)
				if (errno != ENOENT)
					return 0;

				errno = 0;
			}
		}

		/* Read next dir entry in /proc */
		procDirEnt = readdir(procDir);
	}

	closedir(procDir);

	/* Check exit status. */
	if (errno != 0)
		return 0;

	/* Get resources for pid. */
	pidNode = getPidNode(pidHash, pid);

	if (pidNode != NULL) {
		procSnapshot->rss = pidNode->rss;
		procSnapshot->vsize = pidNode->vsize;
		procSnapshot->nproc = 1;

	} else {
		/* No process, exit with 0 counts */
		return 1;
	}

	free(pidNode);

	/* Add child resources for all pids. */
	appendPidQueue(pidQueue, pid);

	while ((pidEntry = takePidQueue(pidQueue)) >= 0) {
		hashIndex = getPidHashIndex(pidHash, pidEntry);

		if (hashIndex >= 0) {
			pidNode = pidHash->hashArray[hashIndex].child;

			while (pidNode != NULL) {
				appendPidQueue(pidQueue, pidNode->pid);

				procSnapshot->rss += pidNode->rss;
				procSnapshot->vsize += pidNode->vsize;
				procSnapshot->nproc += 1;

				pidNode = pidNode->next;
			}
		}
	}

	return 1;
}

/*
 * Determine if a directory name in /proc belongs to a pid. A pid directory must contain only
 * numbers.
 */
int isPidDir(const char *charPtr) {

	if (*charPtr == '\0')
		return 0;

	while (*charPtr != '\0' && isdigit(*charPtr))
		++charPtr;

	if (*charPtr != '\0')
		return 0;

	return 1;
}

/* Function: printSummary */
void printSummary(ProcMonitorSummary *monitorSummary) {

	const char *suffixes[] = {"kB", "MB", "GB", "TB", "EB"};
	const int suffix_size = 5;
	int suffixIndex = 0;

	float vssSize;

	int s;
	int ms;

	if (monitorSummary == NULL)
		return;

	/* Get time. */
	s = monitorSummary->runTime / 1000;
	ms = monitorSummary->runTime % 1000;

	/* Format VSS. */
	vssSize = monitorSummary->maxRss;

	while (vssSize > 1024.0 && suffixIndex < suffix_size) {
		vssSize /= 1024.0;
		++suffixIndex;
	}

	printf("\nRun time: %d.%d (%d:%d:%d.%d)\n", s, ms, s / 3600, (s % 3600) / 60, s % 60, ms);
	printf("Max RSS: %ld kB (%0.2f%s)\n", monitorSummary->maxRss, vssSize, suffixes[suffixIndex]);

	printf("Process exited with code %d\n", monitorSummary->exitStatus);

	return;
}

/* Function: err */
void err(const char *msg, int errCode) {

	fprintf(stderr, "%s: %s\n", progName, msg);

	exit(errCode);
}

/* Function: errnoErr */
void errnoErr(const char *msg, int errCode) {

	fprintf(stderr, "%s: %s: %s\n", progName, msg, strerror(errno));

	exit(errCode);
}

/* Function: printHelp */
void printHelp(void) {
	printf("%s -o <OUT_FILE> [-w <SECS>] <COMMAND>\n", progName);
	printf("%s -h\n\n", progName);

	printf("-h\n");
	printf("\tPrint help information.\n\n");

	printf("-o <OUT_FILE>\n");
	printf("\tOutput CSV file to write.\n\n");

	printf("-w <SECS>\n");
	printf("\tWait <SECS> seconds between each check on the child process being monitored (default = 1).\n\n");
	
	printf("-v\n");
	printf("\tPrint extra (verbose) information about the process after it terminates.\n\n");

	printf("<COMMAND> is the command to be run and monitored with any arguments for the\n");
	printf("command following immediately from it.\n\n");

	return;
}
