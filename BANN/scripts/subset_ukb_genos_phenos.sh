#!/bin/bash

#$ -l h_vmem=32G
#$ -l h_rt=01:00:00
#$ -cwd
#$ -j y

geno_dir=../../data/processed/genotypes
opt=../../opt

source /broad/software/scripts/useuse
reuse -q Anaconda3
source activate ../../opt/for_James/bgen

use R-4.1
R --vanilla <<EOF
library(tidyverse)

set.seed(123)

train_df <- read_csv("../../data/processed/ukb_training_set.csv")
train_subset_df <- train_df %>%
  filter(!is.na(alt_log_adj)) %>%
  select(id, alt_log_adj)
train_subset_df %>%
  write_tsv("../../data/processed/bann/bann_alt_train_phenos.txt")
train_subset_df %>%
  select(fid = id, iid = id) %>%
  write_tsv("../../data/processed/genotypes/bann_train_ids.txt", col_names = FALSE)

test_df <- read_csv("../../data/processed/ukb_testing_set.csv") %>%
  filter(!is.na(alt_log_adj))
test_df %>%
  select(fid = id, iid = id) %>%
  write_tsv("../../data/processed/genotypes/bann_test_ids.txt", col_names = FALSE)
EOF

${opt}/plink2 \
   --pfile ${geno_dir}/ukb_pathway_filtered \
   --keep ${geno_dir}/bann_test_ids.txt \
   --export A \
   --out ${geno_dir}/bann_test_ukb_pathway_filtered
