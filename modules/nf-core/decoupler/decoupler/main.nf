process DECOUPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc0ee6d6033b9f04c6377ad3b1cf5f924e3243626ab8d7be836d9d6617f8da4e/data' :
        'community.wave.seqera.io/library/decoupler-py_matplotlib_pandas_scanpy:369a5afe315b9b30' }"

    input:
    tuple val(meta), path(mat)
    path(net)
    path(gtf)

    output:
    tuple val(meta), path("*estimate_decoupler.tsv"), emit: dc_estimate
    tuple val(meta), path("*pvals_decoupler.tsv"), emit: dc_pvals
    tuple val(meta), path("*decoupler_plot.png"), emit: png
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'decoupler.py'

    stub:
    """
    touch ${task.ext.prefix}_estimate_decoupler.tsv
    touch ${task.ext.prefix}_pvals_decoupler.tsv
    touch ${task.ext.prefix}_decoupler_plot.png
    touch versions.yml
    """
}
