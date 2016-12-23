set -e

# Python 3 (default) environment
. bin/activate kescases;
pip install --upgrade pip;
pip install numpy==1.11.2;
pip install pandas==0.19.1;
pip install pysam==0.9.1.4;
pip install snakemake==3.8.2;
pip install biopython==1.68
. bin/deactivate;

# Python 2 environment
. bin/activate kespy2
pip install pandas==0.19.1;
pip install biopython==1.68
. bin/deactivate kespy2


# Mark complete
touch kescases_pkg_flag
