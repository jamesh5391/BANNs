import numpy as np 
from BANN import *
import os
from utils import rep_col

def loadModel(input_path):
    w = np.load(os.path.join(input_path, 'SNP_weights.npy'))
    pip = np.load(os.path.join(input_path, 'SNP_pips.npy'))
    kernel = np.load(os.path.join(input_path, 'SNP_kernel.npy'))

    return w, pip, kernel

def loadMask(input_path):
    
def main():
    mask = loadMask("../../../data/processed/bann")
    weights, pip, kernel = loadModel("../../../data/processed/bann/model")

    b=rep_col(np.sum(weights * pip * kernel, axis=1),self.mask.shape[1])
	G=np.matmul(self.X,self.mask*b)
if __name__ == "__main__":
    main()
