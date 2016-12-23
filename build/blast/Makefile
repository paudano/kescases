VERSION=2.5.0
BUILD_DIR=ncbi-blast-$(VERSION)+
INSTALL_TARGET=../../bin/blastn
INSTALL_SOURCE=../build/blast/$(BUILD_DIR)/bin
DL_FILE=ncbi-blast-$(VERSION)+-x64-linux.tar.gz
DL_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(VERSION)/ncbi-blast-$(VERSION)+-x64-linux.tar.gz

.SECONDEXPANSION:

$(INSTALL_TARGET): $(BUILD_DIR)/bin/$$(notdir $$@)
	ln -sf $(INSTALL_SOURCE)/$(notdir $@) $@

$(BUILD_DIR)/bin/blastn: $(DL_FILE)
	tar -zxvf $<
	touch $@

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)
