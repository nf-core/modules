import scimap as sm
import anndata as ad
import pandas as pd

filepath = ['$csv_path']
adata = sm.pp.mcmicro_to_scimap(filepath)

#need to save adata
adata.write('scimap_output.h5ad')
