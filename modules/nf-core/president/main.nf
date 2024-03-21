process PRESIDENT {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/president:0.6.8--pyhdfd78af_0' :
        'biocontainers/president:0.6.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path reference
    val compress

    output:
    tuple val(meta), path("output/*.fasta*"), emit: fasta
    tuple val(meta), path("output/*.tsv")   , emit: report
    tuple val(meta), path("output/*.log")   , emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: ''
    if (prefix != ''){
        prefix = "--prefix " + prefix
    }

    """
    mkdir output
    president \\
        --query $fasta \\
        --reference $reference \\
        --path output \\
        --threads $task.cpus \\
        $prefix \\
        $args

    if [ "$compress" = true ] ; then
        gzip output/*.fasta;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        president: \$(president --version |& sed '1!d ; s/president v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: ''
    """
    touch ${prefix}report.tsv
    touch ${prefix}president_logger.log
    touch ${prefix}valid.fasta
    touch ${prefix}invalid.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        president: \$(president --version |& sed '1!d ; s/president v//')
    END_VERSIONS
    """
}
