#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=00:10:00
#$ -cwd
#$ -j y 

source /broad/software/scripts/useuse
reuse -q Anaconda3
source activate /broad/hptmp/jhu/myenv

python <<EOF_PYTHON_CODE

import pandas as pd 


df = pd.read_csv('../data/processed/ukb_testing_set.csv')
ids = pd.read_csv('../data/processed/genotypes/TOY_test_ids.txt', sep='\t', header=None)
print(ids.tail)
print(ids.head)
print(ids.shape)
print(ids)

ids = ids[ids.columns[0]].tolist()
print("ids to list shape", len(ids))

df = df[['id', 'alt_log']]

df = df[df['id'].isin(ids)]
print(df.head)
print(df.shape)
print(df)

df.to_csv('../data/processed/TOY_test_alt_phenos.txt', index=False, sep='\t')


EOF_PYTHON_CODE

