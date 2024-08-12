#!/usr/bin/env python3

# This prevents numba caching from causing issues
# Relevant for Singularity and Apptainer
import os
tmpdir = os.path.join(os.getcwd(), "tmp")
os.environ["NUMBA_CACHE_DIR"] = tmpdir
os.makedirs(tmpdir, exist_ok=True)

import scanpy
import anndata2ri
import rpy2
import rpy2.robjects as ro
import platform
soupx = ro.packages.importr('SoupX')

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

adata = scanpy.read_h5ad("${filtered}")

adata_pp = adata.copy()
scanpy.pp.normalize_per_cell(adata_pp)
scanpy.pp.log1p(adata_pp)

scanpy.pp.pca(adata_pp)
scanpy.pp.neighbors(adata_pp)
scanpy.tl.leiden(adata_pp, key_added="soupx_groups")

soupx_groups = anndata2ri.py2rpy(adata_pp.obs["soupx_groups"])
del adata_pp

sce = anndata2ri.py2rpy(adata)

adata_unfiltered = scanpy.read_h5ad("${unfiltered}")
sce_unfiltered = anndata2ri.py2rpy(adata_unfiltered)

get_counts = ro.r("function(sce) { assay(sce, 'X') }")

data = get_counts(sce)
data_unfiltered = get_counts(sce_unfiltered)

sc = soupx.SoupChannel(data_unfiltered, data, calcSoupProfile = False)

soupProf = ro.r("function(data) { data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data)) }")(data)
sc = soupx.setSoupProfile(sc, soupProf)
sc = soupx.setClusters(sc, soupx_groups)
sc = soupx.autoEstCont(sc, doPlot = False)
out = soupx.adjustCounts(sc, roundToInt = False)

adata.X = anndata2ri.rpy2py(out).T

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": scanpy.__version__,
        "anndata2ri": anndata2ri.__version__,
        "rpy2": rpy2.__version__,
        "SoupX": soupx.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
