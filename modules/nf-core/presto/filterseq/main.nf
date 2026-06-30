process PRESTO_FILTERSEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.9--pyhdfd78af_0':
        'quay.io/biocontainers/presto:0.7.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_quality-pass.fastq.gz"),  emit: reads
    path "*_command_log.txt" , emit: logs
    tuple val("${task.process}"), val('presto'), eval('FilterSeq.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions
    path "*.tab" , emit: log_tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def log_prefix = reads.baseName.replaceFirst(/\.f(ast)?q$/, '')
    """
    FilterSeq.py quality \\
    -s $reads \\
    --outname ${prefix} \\
    --log ${log_prefix}.log \\
    --nproc ${task.cpus} \\
    $args > ${prefix}_command_log.txt

    ParseLog.py -l ${log_prefix}.log $args2 -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def log_prefix = reads.baseName.replaceFirst(/\.f(ast)?q$/, '')
    """
    touch ${prefix}_quality-pass.fastq.gz \\
        ${prefix}_command_log.txt \\
        ${log_prefix}.log \\
        ${log_prefix}_table.tab
    """
}
