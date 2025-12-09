
process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.8--py39h2de1943_0':
        'biocontainers/whatshap:2.8--py39h2de1943_0' }"

    input:
    tuple val(meta), path(vcf) // channel: [ val(meta), path(vcf) ]
    val(include_tsv_output)    // value:   [ true | false ]
    val(include_gtf_output)    // value:   [ true | false ]
    val(inlude_block_output)   // value:   [ true | false ]

    output:
    tuple val(meta), path("${prefix}.tsv"),                                    emit: tsv,   optional: true
    tuple val(meta), path("${prefix}.gtf"),                                    emit: gtf,   optional: true
    tuple val(meta), path("${prefix}.txt"),                                    emit: block, optional: true
    tuple val(meta), path("${prefix}.log"),                                    emit: log
    tuple val("${task.process}"), val('whatshap'), eval("whatshap --version"), emit: versions_whatshap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def output_tsv   = include_tsv_output  ? "--tsv ${prefix}.tsv"              : ''
    def output_gtf   = include_gtf_output  ? "--gtf ${prefix}.gtf"              : ''
    def output_block = inlude_block_output ? "--block-list ${prefix}.txt" : ''
    """
    whatshap stats \\
        $args \\
        $output_tsv \\
        $output_gtf \\
        $output_block \\
        $vcf \\
        | tee ${prefix}.log
    """

    stub:
    def args            = task.ext.args       ?: ''
    prefix              = task.ext.prefix     ?: "${meta.id}"
    def tsv_touch_cmd   = include_tsv_output  ? "touch ${prefix}.tsv" : ''
    def gtf_touch_cmd   = include_gtf_output  ? "touch ${prefix}.gtf" : ''
    def block_touch_cmd = inlude_block_output ? "touch ${prefix}.txt" : ''
    def log_touch_cmd   = "touch ${prefix}.log"
    """
    echo $args

    $tsv_touch_cmd
    $gtf_touch_cmd
    $block_touch_cmd
    $log_touch_cmd
    """
}
