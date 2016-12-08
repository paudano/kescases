VERSION=0.7.15
BUILD_DIR=bwa-$(VERSION)
INSTALL_TARGET=../../bin/bwa
INSTALL_SOURCE=../build/bwa/$(BUILD_DIR)/bwa
DL_FILE=bwa-$(VERSION).tar.bz2
DL_URL=http://downloads.sourceforge.net/project/bio-bwa/$(DL_FILE)

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

