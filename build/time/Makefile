../../bin/time:
	wget https://ftp.gnu.org/gnu/time/time-1.7.tar.gz
	tar -zxf time-1.7.tar.gz
	cd time-1.7; ./configure
	cd time-1.7; make
	mkdir -p $(dir $@)
	ln -sf ../build/time/time-1.7/time $@

.PHONY: clean
clean:
	rm -rf time-1.7*
	rm ../../bin/time

