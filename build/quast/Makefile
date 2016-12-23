VERSION=4.4
BUILD_DIR=quast-$(VERSION)
INSTALL_TARGET=../../bin/quast.py
INSTALL_SOURCE=../build/quast/$(BUILD_DIR)
DL_FILE=quast-$(VERSION).tar.gz
DL_URL=https://downloads.sourceforge.net/project/quast/quast-$(VERSION).tar.gz

.SECONDEXPANSION:

$(INSTALL_TARGET): $(BUILD_DIR)/$$(notdir $$@)
	ln -sf $(INSTALL_SOURCE)/$(notdir $@) $@

$(BUILD_DIR)/quast.py: $(DL_FILE)
	tar -zxvf $<
	touch $@

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

# wget https://downloads.sourceforge.net/project/quast/quast-4.4.tar.gz
