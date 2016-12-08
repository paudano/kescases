REPO_BASE=https://github.com/ncbi

VERSION_NGS=1.3.0
VERSION_NCBI-VDB=2.8.0
VERSION_NGS-TOOLS=2.8.0
VERSION_SRA-TOOLS=2.8.0-1
VERSION_HDF5=1.8.17

DEPS_NGS=
DEPS_NCBI-VDB=jar/ngs-java.jar
DEPS_NGS-TOOLS=lib64/libncbi-vdb.a
DEPS_SRA-TOOLS=bin/sra-search lib/libhdf5.so

CONF_FLAG_NGS=
CONF_FLAG_NCBI-VDB=
CONF_FLAG_NGS-TOOLS=
CONF_FLAG_SRA-TOOLS=--with-hdf5-prefix=$(abspath build/lib)

.SECONDEXPANSION:

# Link
../../bin/fastq-dump: bin/fastq-dump
	ln -s ../build/sra/bin/fastq-dump $@

# Build sra-tools
bin/fastq-dump: build/sra-tools/conf-flag
	cd build/sra-tools; make; make install

# Build ngs-tools
bin/sra-search: build/ngs-tools/conf-flag
	cd build/ngs-tools; make; make install

# Build ncbi-vdb
lib64/libncbi-vdb.a: build/ncbi-vdb/conf-flag
	cd build/ncbi-vdb; make; make install

# Build ngs
jar/ngs-java.jar: build/ngs/conf-flag
	cd build/ngs; make; make install

# Build hdf5
lib/libhdf5.so: build/hdf5/conf-flag
	cd build/hdf5/hdf5-$(VERSION_HDF5); make; make install

### Download ###
build/%/conf-flag: $$(DEPS_$$(shell echo $$* | tr a-z A-Z))
	$(eval NAME=$(shell echo $* | tr a-z A-Z))
	$(eval VERSION=$(VERSION_$(NAME)))
	
	mkdir -p build
	cd build; git clone $(REPO_BASE)/$*.git
	cd build/$*; git checkout $(VERSION)
	cd build/$*; ./configure --prefix=$(abspath .) --build-prefix=$(abspath build/ncbi-outdir) $(CONF_FLAG_$(shell echo $* | tr a-z A-Z))
	touch build/$*/conf-flag

clean:
	rm -rf build

### Download HDF5 Library ###
build/hdf5/conf-flag:
	wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$(VERSION_HDF5)/src/hdf5-$(VERSION_HDF5).tar.gz -P build/hdf5/
	cd build/hdf5/; tar -zxvf hdf5-$(VERSION_HDF5).tar.gz
	cd build/hdf5/hdf5-$(VERSION_HDF5); ./configure --prefix=$(abspath .)
	touch build/hdf5/conf-flag

