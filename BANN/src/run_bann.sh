#!/bin/sh

#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=12:00:00

#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd

bann_dir=../opt/BANN/
X_file=../data/processed/bann/bann_genos.txt
y_file=../data/processed/bann/bann_alt_phenos_ordered.txt
mask_file=../data/processed/bann/bann_mask_ordered.txt

reuse -q Anaconda3

# Navigate to your project directory (where your Python script is)
cd ${bann_dir}/src

python main.py "${X_file}" "${y_file}" "${mask_file}" 