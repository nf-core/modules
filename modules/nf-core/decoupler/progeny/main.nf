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
    tuple val(meta), path("progeny_mlm_pvals_decoupler.tsv")    , emit: dc_pvals
    tuple val(meta), path("acts_decoupler.tsv")                 , emit: dc_acts
    tuple val(meta), path("figures/*umap.png")                  , emit: dc_umap
    tuple val(meta), path("figures/*violin.png")                , emit: dc_violin
    tuple val(meta), path("matrixplot_decoupler.png")           , emit: dc_matrixplot
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'progeny.py'

    stub:
    """
    touch progeny_mlm_estimate_decoupler.tsv
    touch progeny_mlm_pvals_decoupler.tsv
    touch acts_decoupler.tsv
    mkdir -p figures
    touch figures/dc_umap.png
    touch figures/dc_violin.png
    touch matrixplot_decoupler.png
    touch versions.yml
    """

}
