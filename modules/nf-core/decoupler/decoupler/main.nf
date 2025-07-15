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
    tuple val(meta), path("*estimate__decoupler.tsv"), emit: dc_estimate
    tuple val(meta), path("*pvals__decoupler.tsv"), emit: dc_pvals
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'decoupler.py'

    stub:
    """
    touch mlm_estimate__decoupler.tsv
    touch mlm_pvals__decoupler.tsv
    touch versions.yml
    """
}
