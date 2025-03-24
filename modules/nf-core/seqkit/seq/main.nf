process SEQKIT_SEQ {
    tag "$meta.id"
    label 'process_low'
    // File IO can be a bottleneck. See: https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}.*")    , emit: fastx
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def extension   = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    def call_gzip   = extension.endsWith('.gz') ? "| gzip -c $args2" : ''
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    seqkit \\
        seq \\
        --threads $task.cpus \\
        $args \\
        $fastx \\
        $call_gzip \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}"
    def extension   = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
