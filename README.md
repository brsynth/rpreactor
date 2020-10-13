[![Anaconda-Server Badge](https://anaconda.org/brsynth/rpreactor/badges/version.svg)](https://anaconda.org/brsynth/rpreactor)
[![Anaconda-Server Badge](https://anaconda.org/brsynth/rpreactor/badges/license.svg)](https://anaconda.org/brsynth/rpreactor)
[![GitHub Action (flake8 + pytest)](https://github.com/brsynth/rpreactor/workflows/Test%20suite%20(flake8%20+%20pytest)/badge.svg)](https://github.com/brsynth/rpreactor/actions?query=workflow%3A%22Test+suite+%28flake8+%2B+pytest%29%22)

# rpreactor

**A command-line and python package to handle biochemical reaction rules.**

*rpreactor* is designed to use reaction rules from [RetroRules](https://retrorules.org/), 
and to be at the core of more complex bioretrosynthesis tools such as [RetroPathRL](https://github.com/brsynth/RetroPathRL).
It relies extensively on [RDKit](https://rdkit.org/) to handle chemicals and reactions.

Please submit your questions or any new issue you may encounter with *rpreactor* using [GitHub's issue system](https://github.com/brsynth/RetroPathRL/issues).

If you are interested by *rpreactor* you may also be interested in:
* [Reactor](https://chemaxon.com/products/reactor): ChemAxon's "*A high performance virtual synthesis engine*"
* [ATLAS](https://lcsb-databases.epfl.ch/pathways/atlas/): "*A repository of all possible biochemical reactions for synthetic biology and metabolic engineering studies*" by our friends from the [LCSB](https://www.epfl.ch/labs/lcsb/)

## Installation

**Important**: rpreactor needs python version prior to 3.8 and rdkit version 2019.03.1.0.

We strongly recommend you to use [conda package manager](https://docs.conda.io/en/latest/), and to follow those steps: 

```bash
# installation in a new conda environment <myenv>
conda create --name <myenv> -c conda-forge -c brsynth rpreactor
conda activate <myenv>
```

```bash
# installation in an already existing environment <myenv>
conda activate <myenv>
conda install --channel brsynth rpreactor
``` 

## Usage

From command line:
```bash
conda activate <myenv>
python -m rpreactor.cli --help
python -m rpreactor.cli --with_hs true inline --inchi "InChI=1/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)" --rsmarts "([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])"
```

From within a script:
```python
import json
import rpreactor

inchi = 'InChI=1/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
rsmarts = '([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])'

o = rpreactor.RuleBurner(rsmarts_list=[rsmarts], inchi_list=[inchi], with_hs=True)
o.compute()
res = json.loads('[' + ', '.join(o._json) + ']')
```

## For developers

### Development installation

After a git clone:

```bash
cd <repository>
conda env create -f environment.yml -n <dev_env>
conda activate <dev_env>
conda develop -n <dev_env> .
```

You may be prompted to install *conda-build* in your base environment (`conda install conda-build`).
The default conda environment name will be `dev_rpreacor` if not specified by `-n <dev_env>`.

Test your installation with:

```bash
conda activate <dev_env>
python -m rpreactor.cli -h
```

To uninstall:

```bash
conda deactivate
conda env remove -n <dev_env>
```

### Test suite

```bash
cd <repository>
pytest
```

### Build and deployment

The process is automated with GitHub's Action.

If you want to check the build process locally:

```bash
CONDA_BLD_PATH=<repository>/conda-bld
mkdir -p ${CONDA_BLD_PATH} 
cd <repository>

conda env create -f recipe/conda_build_env.yaml -n <build_env>
conda activate <build_env>
conda build -c conda-forge --output-folder ${CONDA_BLD_PATH} recipe

conda convert --platform osx-64 --platform linux-64 --platform win-64 --output-dir ${CONDA_BLD_PATH} ${CONDA_BLD_PATH}/*/rpreactor-*
```

### Develop with Docker

*rpreactor* only works on Linux systems. You may use Docker and [miniconda3 image](https://hub.docker.com/r/continuumio/miniconda3) 
to develop on other operating systems:

```bash
cd <repository>
docker pull continuumio/miniconda3
docker run -d --name dev_rpreactor --mount type=bind,source="$(pwd)",target=/workdir  -it continuumio/miniconda3 bash
```

And then follow the "development installation" steps.
