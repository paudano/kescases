../../bin/traceproc: traceproc
	mkdir -p $(dir $@)
	ln -sf ../build/traceproc/$< $@

traceproc: src/traceproc.c
	gcc -std=gnu99 -lm -lpthread -lrt -o traceproc src/traceproc.c

clean:
	rm traceproc
	rm ../../bin/traceproc

