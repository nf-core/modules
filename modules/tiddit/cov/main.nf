process TIDDIT_COV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::tiddit=2.12.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiddit:2.12.1--py38h1773678_0' :
        'quay.io/biocontainers/tiddit:2.12.1--py38h1773678_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*.tab"), optional: true, emit: cov
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
        --bam $bam \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.wig
    touch ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """
}
