INSTALLER_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

INSTALL_TARGET_DIR=../../bin
INSTALL_TARGET_PREFIX=../build/miniconda/bin
INSTALL_TARGET_ENV_PREFIX=../build/miniconda/envs/kescases/bin
INSTALL_TARGET_ENV2_PREFIX=../build/miniconda/envs/kespy2/bin

ALL_PKGS=$(addprefix envs/kescases/lib/python3.5/site-packages/,pysam pandas numpy) $(addprefix envs/kescases/bin/,snakemake)

.PHONY: all
all: $(INSTALL_TARGET_DIR)/python $(INSTALL_TARGET_DIR)/python2 $(INSTALL_TARGET_DIR)/snakemake


$(INSTALL_TARGET_DIR)/python: kescases_pkg_flag
	ln -sf $(INSTALL_TARGET_ENV_PREFIX)/$(notdir $@) $@

$(INSTALL_TARGET_DIR)/snakemake: kescases_pkg_flag
	ln -sf $(INSTALL_TARGET_ENV_PREFIX)/$(notdir $@) $@

$(INSTALL_TARGET_DIR)/python2: kescases_pkg_flag
	ln -sf $(INSTALL_TARGET_ENV2_PREFIX)/$(notdir $@) $@

kescases_pkg_flag: envs/kescases/bin/python3.5 envs/kespy2/bin/python2.7
	./installpkgs.sh

envs/kescases/bin/python3.5: bin/conda
	bin/conda create --yes --name kescases python=3.5.2
	touch envs/kescases/bin/python3.5

envs/kespy2/bin/python2.7: bin/conda
	bin/conda create --yes --name kespy2 python=2.7.13
	touch envs/kespy2/bin/python2.7

bin/conda: $(notdir $(INSTALLER_URL))
	bash $< -f -b -p `pwd`
	touch bin/conda

$(notdir $(INSTALLER_URL)):
	wget $(INSTALLER_URL)
	touch $@

.PHONY: clean
clean:
	for FILE in $$(/bin/ls); do if [ $$FILE != Makefile -a $$FILE != installpkgs.sh ]; then rm -rf $$FILE; fi; done

