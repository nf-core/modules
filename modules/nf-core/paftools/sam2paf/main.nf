process PAFTOOLS_SAM2PAF {
    tag "$meta.id"
    label 'process_low'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), file("*.paf")          , emit: paf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    samtools view -h ${bam} | paftools.js sam2paf - > ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paftools.js: \$(paftools.js --version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paftools.js: \$(paftools.js --version)
    END VERSIONS
    """
}
