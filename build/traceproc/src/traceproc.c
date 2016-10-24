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


/* Constants. */
const int ERR_NONE = 0;
const int ERR_USAGE = 1;
const int ERR_IO = 2;
const int ERR_THREAD = 3;
const int ERR_PARSE = 4;

const int STAT_BUFFER_SIZE = 1024;

/* Structures */
struct proc_monitor_t {
	long run_time;    /* Run time (ms). */
	long max_rss;     /* Max resident size (memory usage) in kB. */
	int exit_status;  /* Process return code. */
};

typedef struct proc_monitor_t proc_monitor_t;

struct trace_rss_t {
	int pid;         /* Process to monitor. */
	FILE *out_file;  /* File time snapshots are written to. */
	long start_ms;   /* Start time (ms) since epoch. */
	int wait_time;   /* Time to wait (s) between snapshots. */
	int active;      /* Set to FALSE (zero) to tell the trace process to stop. */
	long max_rss;    /* Return max resident size (memory usage) in kB. */
};

typedef struct trace_rss_t trace_rss_t;

/* Function prototypes. */
void monitorProcess(FILE *outFile, char **argv, int waitTime, proc_monitor_t *monitorInfo, int cpuSet);
long getResidentSize(int pid);
void *traceRss(void *traceRssArgs);

void printSummary(proc_monitor_t *monitorInfo);

void printHelp(void);

void err(const char *msg, int errCode);
void errnoErr(const char *msg, int errCode);

/* Program name. */
const int PROG_NAME_SIZE = 25;
char *progName = NULL;

/* Function: main */
int main(int argc, char **argv) {

	/* Declarations. */
	int option;    /* Option character. */
	int inOpts;    /* Flag to stop option processing when false (0). */
	int waitTime;  /* Number of seconds to wait between checking the child process. */
	int verbose;   /* Flag to output verbose information about the process after it terminates. */

	char *outFileName;  /* Name of the output file. */

	char **endPtr;  /* End pointer for strtol. */

	FILE *outFile;

	proc_monitor_t monitorInfo;

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

	/* Set program name for error reporting. */
	if (argc > 0)
		strncpy(progName, basename(argv[0]), PROG_NAME_SIZE);
	else
		strncpy(progName, "trackmem", PROG_NAME_SIZE);

	/* Process arguments. */
	while (inOpts == 1 && (option = getopt(argc, argv, "+c:ho:w:")) != -1) {

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

		fprintf(outFile, "#n\ttime_ms\tmem_kb\n");

	} else {
		outFile = NULL;
	}

	/* Fork and monitor. */
	monitorProcess(outFile, argv + optind, waitTime, &monitorInfo, cpuSet);

	/* Close file. */
	if (outFile != NULL)
		fclose(outFile);

	/* Print summary. */
	if (verbose)
		printSummary(&monitorInfo);

	return monitorInfo.exit_status;
}

/* Function: monitorProcess */
void monitorProcess(FILE *outFile, char **argv, int waitTime, proc_monitor_t *monitorInfo, int cpuSet) {

	/* Declarations. */

	struct timespec currentTime;  /* Timespec struct for getting the current time. */
	long startMs;                 /* Clock time: milliseconds at start. */
	long ms;                      /* Elapsed time: milliseconds. */

	int childPid;  /* Monitored PID. */
	int exitStatus;  /* Monitored PID exit status. */

	trace_rss_t traceArgs;

	pthread_t traceThread;

	cpu_set_t cpuMask;

	int index;

	/* Fork child process. */
	childPid = fork();

	if (childPid == -1)
		errnoErr("Error creating child process to be monitored", ERR_THREAD);

	/* Child: Exec */
	if (childPid == 0) {

		/* Set CPUS */
		if (cpuSet != 0) {

			CPU_ZERO(&cpuMask);

			for (index = 0; index < sizeof(unsigned int) * 8; ++index) {
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

	startMs = currentTime.tv_sec * 1000 + currentTime.tv_nsec / 1000000;

	/* Set traceArgs. */
	traceArgs.pid = childPid;
	traceArgs.out_file = outFile;
	traceArgs.start_ms = startMs;
	traceArgs.wait_time = waitTime;
	traceArgs.active = 1;
	traceArgs.max_rss = 0;

	/* Start trace process. */
	errno = pthread_create(&traceThread, NULL, &traceRss, &traceArgs);

	if (errno != 0)
		errnoErr("Error creating trace thread", ERR_THREAD);

	/* Wait for process to stop. */
	waitpid(childPid, &exitStatus, 0);

	clock_gettime(CLOCK_MONOTONIC, &currentTime);

	ms = currentTime.tv_sec * 1000 + currentTime.tv_nsec / 1000000;

	/* Wait for trace process to stop. */
	traceArgs.active = 0;

	pthread_join(traceThread, NULL);

	/* Set monitor info. */
	if (monitorInfo != NULL) {
		monitorInfo->max_rss = traceArgs.max_rss;
		monitorInfo->run_time = ms - startMs;
		monitorInfo->exit_status = WEXITSTATUS(exitStatus);
	}

	return;
}

/* Function: traceVss */
void *traceRss(void *traceRssArgs) {

	/* Declarations. */
	int count;  /* Number of iterations. */
	long ms;     /* Elapsed time (milliseconds). */

	struct timespec currentTime;  /* Timespec struct for getting the current time. */

	long rss;     /* Resident size (kB) of the process in one cycle. */
	long maxRss;  /* Max RSS. */

	trace_rss_t *traceArgs;

	/* Check arguments. */
	if (traceRssArgs == NULL) {
		pthread_exit(NULL);
	}

	traceArgs = (trace_rss_t *) traceRssArgs;

	/* Monitor (init). */
	ms = traceArgs->start_ms;

	maxRss = 0;
	count = 0;

	clock_gettime(CLOCK_MONOTONIC, &currentTime);
	rss = getResidentSize(traceArgs->pid);

	/* Monitor VSS until the child process dies. */
	while (rss >= 0 && traceArgs->active) {

		/* Set max. */
		if (rss > maxRss)
			maxRss = rss;

		/* Parse time-stamp. */
		ms = currentTime.tv_sec * 1000 + currentTime.tv_nsec / 1000000 - traceArgs->start_ms;

		/* Write. */
		if (traceArgs->out_file != NULL) {
			fprintf(traceArgs->out_file, "%d\t%ld.%ld\t%ld\n", ++count, ms / 1000, ms % 1000, rss);
			fflush(traceArgs->out_file);
		}

		/* Pause between cycles. */
		sleep(traceArgs->wait_time);

		clock_gettime(CLOCK_MONOTONIC, &currentTime);

		rss = getResidentSize(traceArgs->pid);
	}

	traceArgs->max_rss = maxRss;

	/* Exit. */
	pthread_exit(NULL);
}

/* Function: getResidentSize
 *
 * Get the estimated resident size of the process in memory in kB.
 */
long getResidentSize(int pid) {

	/* Declarations. */
	long residentSize;  /* VSS size to be returned. */

	FILE *statFile;  /* Status file (/proc/PID/stat) for the child process. */

	const int STAT_FILE_SIZE = 128;
	char statFileName[STAT_FILE_SIZE];

	const int RET_NO_PID = -1;     /* Return value if the process PID is not alive. */
	const int RET_ERR_ARG = -2;    /* Return value for bad arguments. */
	const int RET_ERR_PARSE = -3;  /* Return value if there is an error parsing stat. */

	char *vssStart;  /* Points to buffer where the VSS field starts. */

	int nRead;  /* Number of bytes read. */

	int offset;

	size_t bufferSize =  STAT_BUFFER_SIZE;
	char *buffer = malloc(sizeof(char) * bufferSize);

	char *sizeStr;

	/* Check arguments. */
	if (pid <= 0) {
		return RET_ERR_ARG;
	}

	/* Open stat. */
	snprintf(statFileName, STAT_FILE_SIZE, "/proc/%d/status", pid);

	statFile = fopen(statFileName, "r");

	if (statFile == NULL) {

		free(buffer);

		return RET_NO_PID;
	}

	/* Get VmRSS line.
	 *
	 * Line format:
	 * VmRSS:      xxxxx kB */
	residentSize = 0;

	while (1) {

		/* Read line. */
		nRead = getline(&buffer, &bufferSize, statFile);

		/* If no more bytes and VmRSS was not found, finish. */
		if (nRead < 0)
			break;

		/* sizeStr will point to the later half (after "VmRSS:") of the line. */
		sizeStr = buffer;

		/* Split at : */
		while (*sizeStr != ':' && sizeStr != '\0')
			++sizeStr;

		if (sizeStr == '\0')
			continue;

		*sizeStr = '\0';

		/* Compare first half to VmRSS. */
		if (strcmp(buffer, "VmRSS") == 0) {

			/* Found the line. sizeStr points to the later
			 * half, use sscanf to convert to a long int. */

			++sizeStr;

			sscanf(sizeStr, "%ld", &residentSize);

			break;
		}
	}

	/* Close and free resources. */
	fclose(statFile);

	free(buffer);

	return residentSize;
}

/* Function: printSummary */
void printSummary(proc_monitor_t *monitorInfo) {

	const char *suffixes[] = {"kB", "MB", "GB", "TB", "EB"};
	const int suffix_size = 5;
	int suffixIndex = 0;

	float vssSize;

	int s;
	int ms;

	if (monitorInfo == NULL)
		return;

	/* Get time. */
	s = monitorInfo->run_time / 1000;
	ms = monitorInfo->run_time % 1000;

	/* Format VSS. */
	vssSize = monitorInfo->max_rss;

	while (vssSize > 1024.0 && suffixIndex < suffix_size) {
		vssSize /= 1024.0;
		++suffixIndex;
	}

	printf("\nRun time: %d.%d (%d:%d:%d.%d)\n", s, ms, s / 3600, (s % 3600) / 60, s % 60, ms);
	printf("Max RSS: %ld kB (%0.2f%s)\n", monitorInfo->max_rss, vssSize, suffixes[suffixIndex]);

	printf("Process exited with code %d\n", monitorInfo->exit_status);

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
