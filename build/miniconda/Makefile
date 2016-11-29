INSTALLER_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

INSTALL_TARGET_DIR=../../bin
INSTALL_TARGET_PREFIX=../build/miniconda/bin
INSTALL_TARGET_ENV_PREFIX=../build/miniconda/envs/kescases/bin

ALL_PKGS=$(addprefix envs/kescases/lib/python3.5/site-packages/,pysam pandas numpy) $(addprefix envs/kescases/bin/,snakemake)

.PHONY: all
all: kescases_pkg_flag

$(INSTALL_TARGET_DIR)/python: build/miniconda/envs/kescases/bin/snakemake
	ln -s $(INSTALL_TARGET_PREFIX)/$(notdir $@) $@

$(INSTALL_TARGET_DIR)/snakemake: 
	ln -s $(INSTALL_TARGET_ENV_PREFIX)/$(notdir $@) $@

kescases_pkg_flag: envs/kescases/bin/python3.5
	./installpkgs.sh

envs/kescases/bin/python3.5: bin/conda
	bin/conda create --yes --name kescases python=3.5.2

bin/conda: $(notdir $(INSTALLER_URL))
	bash Miniconda3-latest-Linux-x86_64.sh -f -b -p `pwd`

$(notdir $(INSTALLER_URL)):
	wget $(INSTALLER_URL)

.PHONY: clean
clean:
	for FILE in $$(/bin/ls); do if [ $$FILE != Makefile -a $$FILE != installpkgs.sh ]; then rm -rf $$FILE; fi; done

