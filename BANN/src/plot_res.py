#!/bin/bash

#$ -l h_vmem=10G
#$ -l h_rt=00:10:00
#$ -cwd
#$ -j y 

source /broad/software/scripts/useuse
reuse -q Anaconda3
source activate /broad/hptmp/jhu/myenv


python <<EOF_PYTHON_CODE

import numpy as np 
import pandas as pd



EOF_PYTHON_CODE
