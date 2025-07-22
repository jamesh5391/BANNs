#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=00:10:00
#$ -cwd
#$ -j y 

source /broad/software/scripts/useuse
reuse -q Anaconda3
source activate /broad/hptmp/jhu/myenv

cd ../opt/BANN/src

python <<EOF_PYTHON_CODE

import numpy as np 
from BANN import *
import os
import pickle 
import pandas as pd
from utils import rep_col


ids_list = "../../../data/processed/genotypes/TOY_test_ids.txt"
pathways_list = "../../../data/processed/prset/TEST_alt_log_prset.snp"


try:
    with open("banns_model.pkl", "rb") as f:
        trained_model = pickle.load(f)
    print("Successfully loaded banns_model.pkl") # Debug print

    X_test = np.loadtxt('../../../data/processed/bann/test/bann_genos.txt')
    print("Successfully loaded X_test.") # Debug print

    mask = np.loadtxt('../../../data/processed/bann/test/bann_mask_ordered.txt')
    print("Successfully loaded mask.") # Debug print

    ids = pd.read_csv(ids_list, header=None, sep='\t').iloc[:,0].tolist()
    pathway_names = pd.read_csv(pathways_list, nrows=0, sep='\t').columns.tolist()
    print("len ids:", len(ids))
    print("len pathways:", len(pathway_names))
    pathway_names = pathway_names[6:]

except FileNotFoundError as e:
    print(f"ERROR: File not found during loading: {e.filename}. Please check your paths.")
    sys.exit(1) # Exit with an error code
except pickle.UnpicklingError as e:
    print(f"ERROR: Could not load banns_model.pkl. File might be corrupted: {e}")
    sys.exit(1)
except Exception as e:
    print(f"ERROR: An unexpected error occurred during initial loading: {e}")
    sys.exit(1)

X_test_scaled = (X_test - trained_model["X_mean"]) / trained_model["X_std"]

#don't worry about biases for now?
print("First mat mul")
SNP_w, SNP_pip, SNP_kernel = trained_model["SNP_w"], trained_model["SNP_pip"], trained_model["SNP_kernel"]
print("SNP_w shape:", SNP_w.shape)
print("SNP_pip shape:", SNP_pip.shape)
print("SNP_kernel shape:", SNP_kernel.shape)

b1=rep_col(np.sum(SNP_w * SNP_pip * SNP_kernel, axis=1), mask.shape[1])
pathway_activation=np.matmul(X_test_scaled, mask*b1)

print("b2 calculation")
SET_w, SET_pip, SET_kernel = trained_model["SET_w"], trained_model["SET_pip"], trained_model["SET_kernel"]
b2=rep_col(np.sum(SET_w * SET_pip * SET_kernel, axis=1), mask.shape[1])

print("y pred")
y_pred_z = np.matmul(pathway_activation, b2)
y_pred = y_pred_z * trained_model["y_std"] + trained_model["y_mean"]


activation_df = pd.DataFrame(pathway_activation, columns=pathway_names, index=ids)
activation_df.to_csv('activations.txt', sep='\t', index=True)
np.savetxt("y_pred.txt", y_pred)

EOF_PYTHON_CODE
