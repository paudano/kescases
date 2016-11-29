VERSION=0.2.6

BUILD_DIR=tabix-$(VERSION)

INSTALL_DIR=../../bin
INSTALL_LINK_DIR=../build/tabix/$(BUILD_DIR)

DL_FILE=tabix-$(VERSION).tar.bz2
DL_URL=http://downloads.sourceforge.net/project/samtools/tabix/tabix-$(VERSION).tar.bz2

TARGETS=tabix bgzip

# Rules

.PHONY: all
all : $(addprefix $(INSTALL_DIR)/, $(TARGETS))

$(INSTALL_DIR)/% : $(BUILD_DIR)/%
	ln -sf $(INSTALL_LINK_DIR)/$* $@

$(addprefix $(BUILD_DIR)/, $(TARGETS)) : $(BUILD_DIR)/Makefile
	cd $(BUILD_DIR); make

$(BUILD_DIR)/Makefile : $(DL_FILE)
	tar -jxvf $<
	touch $@

$(DL_FILE) :
	wget $(DL_URL)

.PHONY: clean
clean :
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

