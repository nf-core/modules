process TRINITY {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.15.2--pl5321hdcf5f25_1':
        'biocontainers/trinity:2.15.2--pl5321hdcf5f25_1' }"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*", arity: '1..*')

    output:
    tuple val(meta), path("*.fa.gz")    , emit: transcript_fasta
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def reads1 = [], reads2 = []
    meta.single_end ? reads1 = reads : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }

    if (meta.single_end) {
        reads_args = "--single ${reads1.join(',')}"
    } else {
        reads_args = "--left ${reads1.join(',')} --right ${reads2.join(',')}"
    }

    // --seqType argument, fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    seqType_args = reads1[0] ==~ /(.*fasta(.gz)?$)|(.*fa(.gz)?$)/ ? "fa" : "fq"

    // Define the memory requirements. Trinity needs this as an option.
    def avail_mem = 7
    if (!task.memory) {
        log.info '[Trinity] Available memory not known - defaulting to 7GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    """
    # Note that Trinity needs the word 'trinity' in the outdir

    Trinity \\
        --seqType ${seqType_args} \\
        --max_memory ${avail_mem}G \\
        ${reads_args} \\
        --output ${prefix}_trinity \\
        --CPU $task.cpus \\
        $args \\
        > >(tee ${prefix}.log)

    gzip \\
        -cf \\
        ${prefix}_trinity.Trinity.fasta \\
        > ${prefix}.fa.gz

    rm ${prefix}_trinity.Trinity.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(Trinity --version | grep 'Trinity version:' | sed 's/Trinity version: Trinity-//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    gzip ${prefix}.fa
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(Trinity --version | grep 'Trinity version:' | sed 's/Trinity version: Trinity-//')
    END_VERSIONS
    """
}
