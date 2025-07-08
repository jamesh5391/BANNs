
import pandas as pd
import numpy as np 

import os

def merge():
    print(os.getcwd())
    sampleDF = pd.read_csv('ukb27892_imp_chrAUT_v3_s487395.sample', sep='\s+', header=0, skiprows=[1])
    sampleDF = sampleDF[['ID_2']]
    sampleDF = sampleDF.rename(columns={'ID_2':'id'})

    print(sampleDF.head(10))

    altDF = pd.read_csv('ukb_validation_set.csv')
    altDF = altDF[['id', 'alt_log']]


    print(altDF.head(10))

    merged = pd.merge(
        sampleDF, 
        altDF,
        on=['id'],
        how='left'
    )

    print(merged.head(100))
    merged = merged.replace([None, '', 'NaN'], np.nan)
    merged = merged.fillna('NA')
    merged.to_csv('mergedALT_ukb_validation_set.txt', sep='\t', index=False)
def initializeX(bgen_path, p):
    #read in bgen 
    #read in phenotype dataset and get all IDs

