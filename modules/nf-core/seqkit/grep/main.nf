process SEQKIT_GREP {
    tag "$meta.id"
    label 'process_low'


    conda "bioconda::seqkit=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0':
        'biocontainers/seqkit:2.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(pattern)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.{fa,fq}.gz")  , emit: fasta
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def suffix = task.ext.suffix ?: "${reference}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"

    """
    seqkit \\
        grep \\
        $args \\
        --threads $task.cpus \\
        -f ${pattern} \\
        ${reference} \\
        -o ${prefix}.${suffix}.gz \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
