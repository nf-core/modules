process DECOUPLER_PROGENY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d304d1b25aa80ce8c44c8d34ea45d6b4f6f50697a4effcf1b9be4a54db19928/data' :
        'community.wave.seqera.io/library/decoupler-py_matplotlib_pandas_pip_pruned:0d4681dad9987ec5' }"

    input:
    tuple val(meta), path(mat)
    val(organism)
    val(top_genes)

    output:
    tuple val(meta), path("progeny_mlm_estimate_decoupler.tsv") , emit: dc_estimate
    tuple val(meta), path("progeny_mlm_pvals_decoupler.tsv"), emit: dc_pvals
    tuple val(meta), path("acts_decoupler.tsv"), emit: dc_acts
    tuple val(meta), path("figures/*umap.png"), emit: dc_umap
    tuple val(meta), path("figures/*violin.png"), emit: dc_violin
    tuple val(meta), path("matrixplot_decoupler.png"), emit: dc_matrixplot
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "{}"
"""
    #!/usr/bin/env python3
    import os

    os.environ["OMNIPATH_CACHE_DIR"] = "./tmp"
    os.makedirs(os.environ["OMNIPATH_CACHE_DIR"], exist_ok=True)

    os.environ["MPLCONFIGDIR"] = "./tmp"
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import pandas as pd
    import matplotlib.pyplot as plt
    import scanpy as sc
    import decoupler as dc

    adata = sc.read_h5ad('${mat}')

    progeny = dc.get_progeny(organism='${organism}', top=${top_genes})

    dc.run_mlm(
        mat=adata,
        net=progeny,
        source='source',
        target='target',
        weight='weight',
        verbose=True
    )

    adata.obsm['progeny_mlm_estimate'] = adata.obsm['mlm_estimate'].copy()
    adata.obsm['progeny_mlm_pvals'] = adata.obsm['mlm_pvals'].copy()

    mlm_estimates = pd.DataFrame(adata.obsm['mlm_estimate'], index=adata.obs_names)
    mlm_pvals = pd.DataFrame(adata.obsm['mlm_pvals'], index=adata.obs_names)
    mlm_estimates.to_csv("progeny_mlm_estimate_decoupler.tsv", sep="\t")
    mlm_pvals.to_csv("progeny_mlm_pvals_decoupler.tsv", sep="\t")

    acts = dc.get_acts(adata, obsm_key='mlm_estimate')

    sc.pl.umap(
        acts,
        color=['Trail', 'louvain'],
        cmap='RdBu_r',
        vcenter=0,
        save="Trail_louvain_umap.png"
    )

    sc.pl.violin(
        acts,
        keys=['Trail'],
        groupby='louvain',
        rotation=90,
        save="Trail_violin.png"
    )

    sc.pl.matrixplot(
        acts,
        var_names=acts.var_names,
        groupby='louvain',
        dendrogram=True,
        standard_scale='var',
        colorbar_title='Z-scaled scores',
        cmap='RdBu_r'
    )
    plt.savefig("matrixplot_decoupler.png", dpi=300)

    acts_df = pd.DataFrame(acts.X, index=acts.obs_names, columns=acts.var_names)
    acts_df.to_csv("acts_decoupler.tsv", sep="\t")

    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
"""

}
