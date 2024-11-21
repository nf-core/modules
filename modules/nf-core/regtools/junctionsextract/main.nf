process REGTOOLS_JUNCTIONSEXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::regtools=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/regtools:1.0.0--452281951f764f78' :
        'community.wave.seqera.io/library/regtools:1.0.0--461ddf16709a70cf' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.junc"), emit: junc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    regtools junctions extract \\
        $args \\
        -s XS \\
        -o ${prefix}.junc \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | grep "Version:" | sed 's/Version:\t//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | grep "Version:" | sed 's/Version:\t//')
    END_VERSIONS
    """
}
