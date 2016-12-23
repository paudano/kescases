VERSION=3.3.2
BUILD_DIR=R-$(VERSION)
INSTALL_TARGET=../../bin/R ../../bin/Rscript
DL_FILE=R-$(VERSION).tar.gz
DL_URL=https://cran.r-project.org/src/base/R-3/R-$(VERSION).tar.gz

.SECONDEXPANSION:

.PHONY: all
all: $(INSTALL_TARGET)

# Install
$(INSTALL_TARGET): bin/$$(notdir $$@) lib64/R/library/ggplot2/R/ggplot2
	ln -sf ../build/r/bin/$(notdir $@) $@

# Install packages: ggplot2
lib64/R/library/ggplot2/R/ggplot2: lib64/R/etc/Rprofile.site
	bin/R -q -e 'install.packages("ggplot2", repos="http://cran.us.r-project.org")'

# Setup site profile
lib64/R/etc/Rprofile.site: bin/R
	echo ".libPaths('$(realpath lib64/R/library)')" > $@

# Install to this directory
bin/R bin/Rscript: $(BUILD_DIR)/bin/$$(notdir $$@)
	cd $(BUILD_DIR); make install

# Build
$(BUILD_DIR)/bin/R $(BUILD_DIR)/bin/Rscript: $(BUILD_DIR)/Makefile
	cd $(BUILD_DIR); make

# Run configure (create Makefile)
$(BUILD_DIR)/Makefile: $(BUILD_DIR)/configure
	cd $(BUILD_DIR); ./configure --prefix=$(realpath .) --with-x=no

# Expand download
$(BUILD_DIR)/configure: $(DL_FILE)
	tar -zxvf $<
	touch $@

# Download file
$(DL_FILE):
	wget $(DL_URL)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(DL_FILE) $(INSTALL_TARGET)

