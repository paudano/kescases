CFLAGS=-Iinclude -std=gnu99 -Iinclude -ggdb
LDLIBS=-lm -lpthread -lrt

../../bin/traceproc : traceproc
	mkdir -p $(dir $@)
	ln -sf ../build/traceproc/$< $@

traceproc : traceproc.o pidhash.o pidqueue.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@
	

%.o : src/%.c include/%.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

traceproc.o : src/traceproc.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean :
	rm traceproc
	rm *.o
	rm ../../bin/traceproc

