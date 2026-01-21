process NGSCHECKMATE_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.1--py312pl5321h577a1d6_4':
        'biocontainers/ngscheckmate:1.0.1--py312pl5321h577a1d6_4' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(snp_pt)

    output:
    tuple val(meta), path("*.vaf"), emit: vaf
    tuple val("${task.process}"), val('ngscheckmate'), eval("ncm.py --help | sed '7!d;s/.* v//g'"), topic: versions, emit: versions_ngscheckmate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq2command  = ( reads instanceof List && reads.size() == 2 ) ? " -2 ${reads[1]} " : ""

    """
    ngscheckmate_fastq -1 ${reads[0]} $fastq2command ${snp_pt} -p ${task.cpus} $args > ${prefix}.vaf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vaf
    """
}
