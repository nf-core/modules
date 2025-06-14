import random

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

# Set random seed
random.seed(52)

# Generate signal data
signal = stats.poisson.rvs(1000, size=990)
doublet_signal = stats.poisson.rvs(1000, size=10)

# Generate background data
x = np.reshape(stats.poisson.rvs(500, size=10000), (1000, 10))

# Insert signal
for idx, signal_count in enumerate(signal):
    col_pos = idx % 10
    x[idx, col_pos] = signal_count

# Insert doublet signal
for idx, signal_count in enumerate(doublet_signal):
    col_pos = (idx % 10) - 1
    x[idx, col_pos] = signal_count

# Create AnnData object
obs_df = pd.DataFrame(x, columns=[str(i) for i in range(x.shape[1])])
adata = AnnData(X=np.random.randint(0, 100, size=x.shape), obs=obs_df)
adata.write("hashsolo_anndata.h5ad")
