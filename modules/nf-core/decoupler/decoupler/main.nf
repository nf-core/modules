process DECOUPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/5369022d8f91af004bd03759a30d463f48ba6e915b734938355018069b53e4e4/data' :
        'community.wave.seqera.io/library/decoupler-py_matplotlib_pandas_scanpy_pruned:c1ac1f1e74a97858' }"

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
