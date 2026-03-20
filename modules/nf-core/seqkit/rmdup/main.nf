process SEQKIT_RMDUP {
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
    tuple val(meta), path("${prefix}.${extension}") , emit: fastx
    tuple val(meta), path("*.log")                  , emit: log
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    extension       = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    // SeqKit/rmdup takes care of compressing the output: https://bioinf.shenwei.me/seqkit/usage/#rmdup
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    seqkit \\
        rmdup \\
        --threads $task.cpus \\
        $args \\
        $fastx \\
        -o ${prefix}.${extension} \\
        2>| >(tee ${prefix}.log >&2)
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}"
    extension       = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}
    echo \\
        '[INFO] 0 duplicated records removed' \\
        > ${prefix}.log
    """
}
