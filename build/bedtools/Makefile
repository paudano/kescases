VERSION=2.25.0
BUILD_DIR=bedtools2
INSTALL_TARGET=../../bin
INSTALL_SOURCE=../build/bedtools/$(BUILD_DIR)/bin
DL_FILE=bedtools-$(VERSION).tar.gz
DL_URL=https://github.com/arq5x/bedtools2/releases/download/v2.25.0/$(DL_FILE)

.SECONDEXPANSION:

$(BUILD_DIR)/bin/bedtools: $(BUILD_DIR)/Makefile
	make -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: $(DL_FILE)
	tar -zxvf $<
	touch $@

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)
