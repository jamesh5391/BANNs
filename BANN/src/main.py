import numpy as np
from BANN import *
import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("X_path", type=str)
    parser.add_argument("y_path", type=str)
    parser.add_argument("mask_path", type=str)
    
    args = parser.parse_args()

    X = np.loadtxt(args.X_path)
    y = np.loadtxt(args.y_path)
    mask = np.loadtxt(args.mask_path)

    bann=BANNs(X,y, mask, nModelsSNP=20, nModelsSET=20)
    [SNP_layer, SET_layer]=bann.run()

    SNPpips=SNP_layer.pip
    SETpips=SET_layer.pip

    plt.figure()
    plt.scatter(np.arange(len(SNPpips)), SNPpips)
    plt.savefig("SNPpips.png")
    
    plt.figure()
    plt.scatter(np.arange(len(SETpips)), SETpips)
    plt.savefig("SETpips.png")

    #Export model parameters
    bann.save_model()


    
    
