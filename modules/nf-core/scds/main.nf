process SCDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-scds_r-seurat:9433ff91d563da96':
        'community.wave.seqera.io/library/bioconductor-scds_r-seurat:0a9b33d0ff3ee819' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.csv"), emit: predictions
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'scds.R'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds
    touch ${prefix}.csv
    touch versions.yml
    """
}
