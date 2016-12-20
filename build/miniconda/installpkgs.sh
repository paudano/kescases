set -e

. bin/activate kescases;
pip install --upgrade pip;
pip install numpy==1.11.2;
pip install pandas==0.19.1;
pip install pysam==0.9.1.4;
pip install snakemake==3.8.2;
. bin/deactivate;
touch kescases_pkg_flag

