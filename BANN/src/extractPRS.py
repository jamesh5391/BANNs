import numpy as np 
from BANN import *
import os
import pickle 
from utils import rep_col

def main():
    #Change to load actual test dataset with correct direcotyr path
    X_test = np.loadtxt('../../../data/processed/genotypes/bann_genos_test.txt')
    mask = np.loadtxt('bann_mask_ordered.txt')
    
    with open("banns_model.pkl", "rb") as f: 
        trained_model = pickle.load(f)

    
    X_test_scaled = (X_test - trained_model["X_mean"]) / trained_model["X_std"]

    #don't worry about biases for now?
    SNP_w, SNP_pip, SNP_kernel = trained_model["SNP_w"], trained_model["SNP_pip"], trained_model["SNP_kernel"]
    b1=rep_col(np.sum(SNP_w * SNP_pip * SNP_kernel, axis=1), mask.shape[1])
    G=np.matmul(X_test_scaled, mask*b1)

    SET_w, SET_pip, SET_kernel = trained_model["SET_w"], trained_model["SET_pip"], trained_model["SET_kernel"]
    b2=rep_col(np.sum(SET_w * SET_pip * SET_kernel, axis=1), mask.shape[1])

    y_pred_z = np.matmul(G, b2)
    y_pred = y_pred_z * trained_model["y_std"] + trained_model["y_mean"]

    np.savetxt("activations.txt", G)
    np.savetxt("y_pred.txt", y_pred)