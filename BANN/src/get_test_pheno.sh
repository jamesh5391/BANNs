#!/bin/bash

#$ -l h_vmem=30G
#$ -l h_rt=1:00:00

#$ -cwd
#$ -j y


pheno_dir=../data/processed

ukb_bgen_dir=/broad/ukbb/imputed_v3 
ukb_sample_dir=/humgen/florezlab/UKBB_app27892 


source /broad/software/scripts/useuse
use R-4.1


R --vanilla <<EOF
library(tidyverse)

set.seed(123)

train_df <- read_csv("../data/processed/ukb_training_set.csv")
train_subset_idx <- sample(seq(1, nrow(train_df)), size = 5000, replace = FALSE)
train_subset_df <- train_df %>%
  slice(train_subset_idx) %>%
  filter(!is.na(alt_log)) %>%
  select(id, alt_log)
train_subset_df %>%
  write_tsv("../data/processed/bann/bann_alt_phenos.txt")
train_subset_df %>%
  select(fid = id, iid = id) %>%
  write_tsv("../data/processed/genotypes/train_bann_ids.txt", col_names = FALSE)

test_df <- read_csv("../data/processed/ukb_testing_set.csv") %>%
  filter(!is.na(alt_log))

test_df %>%
  select(fid = id, iid = id) %>%
  write_tsv("../data/processed/genotypes/final_test_ids.txt", col_names = FALSE)
toy_test_subset_idx <- sample(seq(1, nrow(test_df)), size = 1000, replace = FALSE)
test_df %>%
  slice(toy_test_subset_idx) %>%
  select(fid = id, iid = id) %>%
  write_tsv("../data/processed/genotypes/TOY_test_ids.txt", col_names = FALSE)
EOF
