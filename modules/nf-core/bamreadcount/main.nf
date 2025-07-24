process BAMREADCOUNT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bam-readcount:1.0.1--h9aeec6d_3':
        'biocontainers/bam-readcount:1.0.1--h9aeec6d_3' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(bed)

    output:
    tuple val(meta), path("*.rc"), emit: rc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam-readcount \\
        $args \\
        -w 1 \\
        -f ${reference} \\
        -l ${bed} \\
        $bam \\
        > ${prefix}.rc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamreadcount: \$(bam-readcount --version | awk -F'version: ' '{print \$2}' | awk -F'-' '{print \$1}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.rc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamreadcount: \$(bam-readcount --version | awk -F'version: ' '{print \$2}' | awk -F'-' '{print \$1}')
    END_VERSIONS
    """
}
