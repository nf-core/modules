process REGTOOLS_JUNCTIONSEXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::regtools=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6c0653189b95b22e16038f61ade205a865857f54eeae9ba0184490a1834f7c9/data' :
        'community.wave.seqera.io/library/regtools:0.5.0--b9a260c4c898354a' }"

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
        -o ${prefix}.junc \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | sed 's/^.*regtools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | sed 's/^.*regtools //; s/ .*\$//')
    END_VERSIONS
    """
}
