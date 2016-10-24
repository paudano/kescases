.PHONY: all
all: bin/time

bin/time:
	make -C build/time

bin/traceproc:
	make -C build/traceproc

.PHONY: clean
clean:
	make -C build/time clean
	make -C build/traceproc clean

