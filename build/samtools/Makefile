VERSION=1.3.1
BUILD_DIR=samtools-$(VERSION)
INSTALL_TARGET=../../bin/samtools
INSTALL_SOURCE=../build/samtools/$(BUILD_DIR)/samtools
DL_FILE=samtools-$(VERSION).tar.bz2
DL_URL=https://sourceforge.net/projects/samtools/files/samtools/$(VERSION)/$(DL_FILE)

$(INSTALL_TARGET): $(INSTALL_SOURCE)
	ln -sf $< $@

$(INSTALL_SOURCE): $(BUILD_DIR)/Makefile
	cd $(BUILD_DIR); make

$(BUILD_DIR)/Makefile: $(DL_FILE)
	tar -jxvf $<
	touch $@

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

