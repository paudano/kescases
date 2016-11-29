VERSION=3.7.1
BUILD_DIR=rtg-tools-$(VERSION)
INSTALL_TARGET=../../bin/rtg
INSTALL_SOURCE=../build/rtg/$(BUILD_DIR)/rtg
DL_FILE=rtg-tools-$(VERSION)-linux-x64.zip
DL_URL=https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.7.1/rtg-tools-3.7.1-linux-x64.zip

$(INSTALL_TARGET): $(BUILD_DIR)/rtg
	ln -sf $(INSTALL_SOURCE) $@

$(BUILD_DIR)/rtg: $(DL_FILE)
	unzip $(DL_FILE)
	touch $@
	cp rtg.cfg $(BUILD_DIR)

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

