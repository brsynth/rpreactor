# Workflow for reaction rule firing

Thomas Duigou (thomas.duigou@inra.fr), INRA, 2018

## Installation
```
# Build conda env
conda create --name pyrule python=3.6
source activate pyrule
conda install --channel rdkit rdkit=2019.03.1.0
conda install --channel bioconda --channel conda-forge snakemake=5.4.0
```

## Specific installation for migale cluster
```
# Download & install locally conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  # Keep default install path 

# Append to .bashrc to enable conda commands
echo "# Prevent issue with buggy PYTHONPATH" >> ~/.bashrc
echo "unset PYTHONPATH" >> ~/.bashrc
echo "# Enable access to Miniconda3 local install" >> ~/.bashrc
echo "source $HOME/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc

# Follow standard installation above
# ...

# Append specific job handler package
source activate pyrule
conda install drmaa
```

## Usage