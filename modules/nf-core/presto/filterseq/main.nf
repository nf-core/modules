process PRESTO_FILTERSEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0':
        'biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_quality-pass.fastq"),  emit: reads
    path "*_command_log.txt" , emit: logs
    path "versions.yml" , emit: versions
    path "*.tab" , emit: log_tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    FilterSeq.py quality \\
    -s $reads \\
    --outname ${prefix} \\
    --log ${reads.baseName}.log \\
    --nproc ${task.cpus} \\
    $args > ${prefix}_command_log.txt

    ParseLog.py -l ${reads.baseName}.log $args2 -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_quality-pass.fastq ${prefix}_quality-pass.fastq \\
        ${prefix}_command_log.txt ${prefix}.log ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
