VERSION=3.9.0
BUILD_DIR=SPAdes-$(VERSION)
INSTALL_TARGET=../../bin/spades.py
INSTALL_SOURCE=../build/spades/$(BUILD_DIR)/spades.py
DL_FILE=SPAdes-$(VERSION).tar.gz
DL_URL=http://spades.bioinf.spbau.ru/release$(VERSION)/SPAdes-$(VERSION).tar.gz

$(INSTALL_TARGET): $(INSTALL_SOURCE)
	ln -sf $< $@

$(INSTALL_SOURCE): $(BUILD_DIR)/spades_compile.sh
	cd $(BUILD_DIR); ./spades_compile.sh

$(BUILD_DIR)/spades_compile.sh: $(DL_FILE)
	tar -xvf $<
	touch $@

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

