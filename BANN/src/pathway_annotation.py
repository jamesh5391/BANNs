import pandas as pd
import numpy as np 
from tqdm import tqdm #just for progress bars 
import natsort as ns #natural sorting
import sys
import time 
import math

# map snps to genes to pathways 
#original code only goes from snps to gene (which are the snp sets)
#to do: change annotation matrix s.t. it is n x p (where p is num pathways) rather than number of gene sets 


	

'''
INPUTS:
    -GWAS summary statistics to collect all SNPs and their ids 
    -0
'''