process TIDDIT_COV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiddit:3.9.3--py39hff726c5_1' :
        'biocontainers/tiddit:3.9.3--py39hff726c5_1' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bed"), optional: true, emit: cov
    tuple val(meta), path("*.wig"), optional: true, emit: wig
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref $fasta" : ""
    """
    tiddit \\
        --cov \\
        -o $prefix \\
        $args \\
        --bam $input \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.wig
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """
}
