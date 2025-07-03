import pandas as pd
import os
from BANN.src.annotation_snp_gene_path import * 


msigdb_path = "BANN/data/msigdb.v2025.1.Hs.symbols.gmt"

df = read_pathway_file(msigdb_path)

print(df.iloc[0])