process EIDO_VALIDATE {
    tag "$samplesheet"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f89ca27f1ccaa40dfcf8d9f6e6aab5f9599c5e1e37cf694c4e4f4ba0641577d8/data' :
        'community.wave.seqera.io/library/eido_peppy:3721c3f85cc3d076' }"

    input:
    path samplesheet
    path schema

    output:
    path "*.log", emit: log
    tuple val("${task.process}"), val("eido"), eval("eido --version 2>&1 | sed 's/^.*eido //;s/ .*//'"), topic: versions, emit: versions_eido

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "validation"
    """
    eido validate $args $samplesheet -s $schema -e > ${prefix}.log
    """

    stub:
    def prefix = task.ext.prefix ?: "validation"

    """
    touch ${prefix}.log
    """
}
