process SEQKIT_GREP {
    tag "$meta.id"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(sequence)
    path pattern

    output:
    tuple val(meta), path("*.{fa,fq,fa.gz,fq.gz}"), emit: filter
    tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def compression_suffix = sequence.getExtension() == "gz" ? ".gz" : ""
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def output_type = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"
    def suffix = output_type + compression_suffix
    def pattern_file = pattern ? "-f ${pattern}" : ""
    """
    seqkit \\
        grep \\
        $args \\
        --threads $task.cpus \\
        ${pattern_file} \\
        ${sequence} \\
        -o ${prefix}.${suffix} \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def compression_suffix = sequence.getExtension() == "gz" ? ".gz" : ""
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def output_type = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"
    def suffix = output_type + compression_suffix
    """
    echo ${args}
    echo "" | gzip > ${prefix}.${suffix}.gz
    """
}
