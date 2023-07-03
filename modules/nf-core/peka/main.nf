process PEKA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::peka=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/peka:1.0.0--pyhdfd78af_0':
        'biocontainers/peka:1.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(peaks)
    tuple val(meta), path(crosslinks)
    path fasta
    path fai
    path gtf

    output:
    tuple val(meta), path("*mer_cluster_distribution*"), emit: cluster,      optional: true
    tuple val(meta), path("*mer_distribution*")        , emit: distribution, optional: true
    tuple val(meta), path("*.pdf")                     , emit: pdf,          optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def VERSION = '1.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # If the modification date and time of the fai is before the fasta then
    # there will be an error. Touching the file first avoids that.
    touch $fai
    mkdir tmp
    TMPDIR=\$(pwd)/tmp peka \
        -i $peaks \
        -x $crosslinks \
        -g $fasta \
        -gi $fai \
        -r $gtf \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peka: $VERSION
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_4mer_cluster_distribution_genome.tsv
    touch ${prefix}_4mer_distribution_genome.tsv
    touch ${prefix}_4mer_genome.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peka: $VERSION
    END_VERSIONS
    """
}
