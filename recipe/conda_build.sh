PKG_NAME=vast2
USER=tfursten

OS=linux-64
#mkdir ~/conda-bld
conda config --set anaconda_upload no
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
export CONDA_BLD_PATH=~/conda-bld
export VERSION="2.0.0"
conda build ..
conda install vast2 --use-local
