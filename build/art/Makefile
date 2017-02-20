BUILD_DIR=art_bin_MountRainier
INSTALL_TARGET=../../bin/art_illumina
INSTALL_SOURCE=../build/art/$(BUILD_DIR)/art_illumina
DL_FILE=artbinmountrainier20160605linux64tgz.tgz
DL_URL=https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz


$(INSTALL_TARGET): $(BUILD_DIR)/art_illumina
	ln -sf ../build/art/$(BUILD_DIR)/art_illumina ../../bin/art_illumina

$(BUILD_DIR)/art_illumina: $(DL_FILE)
	tar -zxvf $<
	touch $@

$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

