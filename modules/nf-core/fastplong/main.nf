process FASTPLONG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastplong:0.3.0--h224cc79_0':
        'biocontainers/fastplong:0.3.0--h224cc79_0' }"

    input:
    tuple val(meta), path(reads)
    path  adapter_fasta
    val   discard_trimmed_pass
    val   save_trimmed_fail

    output:
    tuple val(meta), path('*.fastplong.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    def fail_fastq  = save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ''
    def output_file = discard_trimmed_pass ? '' : "--out ${prefix}.fastplong.fastq.gz"
    def report_title = task.ext.report_title ?: "${prefix}_fastplong_report"
    """
    fastplong \\
        --in ${reads} \\
        $output_file \\
        --json ${prefix}.fastplong.json \\
        --html ${prefix}.fastplong.html \\
        $adapter_list \\
        $fail_fastq \\
        --thread $task.cpus \\
        --report_title $report_title\\
        $args \\
        2> >(tee ${prefix}.fastplong.log >&2)


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e "s/fastplong //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_reads = discard_trimmed_pass ? "" : "echo '' | gzip > ${prefix}.fastplong.fastq.gz"
    def touch_fail  = save_trimmed_fail ? "echo '' | gzip > ${prefix}.fail.fastq.gz" : ""
    """
    echo $args

    $touch_reads
    $touch_fail
    touch ${prefix}.fastplong.json
    touch ${prefix}.fastplong.html
    touch ${prefix}.fastplong.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e "s/fastplong //g")
    END_VERSIONS
    """
}
