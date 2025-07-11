
import pandas as pd
import numpy as np 
import bgen_reader
import os
import sys

def merge(pheno_filepath):
    print(os.getcwd())
    sampleDF = pd.read_csv('ukb27892_imp_chrAUT_v3_s487395.sample', sep='\s+', header=0, skiprows=[1])
    sampleDF = sampleDF[['ID_2']]
    sampleDF = sampleDF.rename(columns={'ID_2':'id'})

    print(sampleDF.head(10))

    altDF = pd.read_csv(pheno_filepath)
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
    #merged.to_csv('mergedALT_ukb_validation_set.txt', sep='\t', index=False)
    return merged

def initializeXY(bgen_path, sample_path, merged_pheno):
    #get indices of all individuals without missing phenotype values
    n = len(merged_pheno)
    print(f"{n} values in merged phenotype dataset")

    non_missing_pheno_mask = merged_pheno['alt_log'].notna()
    print(f"{n - np.sum(non_missing_pheno_mask)} missing values in merged phenotype dataset")
    indices_to_load = np.where(non_missing_pheno_mask)[0].tolist()

    #read in bgen
    try:
        bgen = bgen_reader.open_bgen(bgen_path, sample_path)
    #filter out all individuals in X with missing phenotype
        X = bgen.read(
            samples_by_idx=indices_to_load,
            dtype=np.float32,
            dosage=True
        
        )
        
    except FileNotFoundError:
        sys.exit(f"Error: BGEN file or sample file not found. Check paths: {bgen_path}, {sample_path}")
    except Exception as e:
        sys.exit(f"Error reading BGEN for specific individuals: {e}")

    #filter out all individuals in X and y with a missing genotype
    non_missing_geno_mask = ~np.isnan(X).any(axis=1)
    X_filtered = X[non_missing_geno_mask, :]

    pheno_mask = non_missing_geno_mask & non_missing_pheno_mask
    y_filtered = merged_pheno['alt_log'][pheno_mask].values()

    assert X_filtered.shape[0] == y_filtered.shape[0], "Number of individuals don't match between X and y!"
    


    return X_filtered, y_filtered


#IMPORTANT!!!! MAKE SURE SNP ORDER OF FILE IS SAME AS X
def initializePathwayMask(annotation_file_path):
    df = pd.read_csv(annotation_file_path, sep='\t')
    df = df.iloc[:, 6:]
    return df.to_numpy()


merged_pheno = merge('ukb_testing_set.csv')
X, y = initializeXY(
    bgen_path = "test.bgen",
    sample_path="sample.test",
    merged_pheno=merged_pheno )
mask = initializePathwayMask("annotation_path.test")


# bann=BANNs(X,y, mask, nModelsSNP=20, nModelsSET=20)
# [SNP_layer, SET_layer]=bann.run()
# print("PVE")
# print(SNP_layer.pve)
# print(SET_layer.pve)


# pips=SNP_layer.pip
# pips2=SET_layer.pip