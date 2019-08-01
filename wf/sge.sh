#$ -S /bin/bash

export PATH=$HOME/miniconda3/bin:$PATH
unset PYTHONPATH
source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/mcts
{exec_job}