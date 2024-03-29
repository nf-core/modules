process VELOCYTO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/velocyto.py:0.17.17--py38h24c8ff8_6':
        'biocontainers/velocyto.py:0.17.17--py38h24c8ff8_6' }"

    stageInMode 'copy'

    input:
    tuple val(meta), path(barcodes), path(bam), path(sorted_bam)
    path gtf

    output:
    tuple val(meta), path("*.loom"), path("*.velocyto.log"), emit: loom
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    velocyto run $args -e ${meta.id} -b ${barcodes} -o . ${bam} ${gtf} > ${prefix}.velocyto.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        velocyto: \$(echo \$(velocyto --version) | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.loom
    touch ${prefix}.velocyto.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        velocyto: \$(echo \$(velocyto --version) | sed 's/^.*version //')
    END_VERSIONS
    """
}
