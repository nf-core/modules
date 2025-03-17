process DECOUPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d304d1b25aa80ce8c44c8d34ea45d6b4f6f50697a4effcf1b9be4a54db19928/data' :
        'community.wave.seqera.io/library/decoupler-py_matplotlib_pandas_pip_pruned:0d4681dad9987ec5' }"

    input:
    tuple val(meta), path(mat)
    path(net)

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
