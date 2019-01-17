#$ -S /bin/bash

unset PYTHONPATH
source $HOME/miniconda3/bin/activate $HOME/.conda/envs/pyrule
{exec_job}
