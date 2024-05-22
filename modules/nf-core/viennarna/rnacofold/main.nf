process VIENNARNA_RNACOFOLD {
    tag "$rnacofold_fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/viennarna:2.6.4--py310pl5321h6cc9453_1':
        'biocontainers/viennarna:2.6.4--py310pl5321h6cc9453_1' }"

    input:
    tuple val(meta), path(rnacofold_fasta)

    output:
    tuple val(meta), path("*.csv")  , emit: rnacofold_csv
    tuple val(meta), path("*.ps")   , emit: rnacofold_ps
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNAcofold < ${rnacofold_fasta} \\
        --jobs=${task.cpus} \\
        --output-format="D" \\
        > ${rnacofold_fasta.baseName}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNAcofold: \$( RNAcofold --version |& sed '1!d ; s/RNAcofold //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${rnacofold_fasta.baseName}.csv
    touch ${rnacofold_fasta.baseName}.ps

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNAcofold: \$( RNAcofold --version |& sed '1!d ; s/RNAcofold //')
    END_VERSIONS
    """
}
