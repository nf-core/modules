process PMDTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmdtools:0.60--hdfd78af_5' :
        'quay.io/biocontainers/pmdtools:0.60--hdfd78af_5' }"

    input:
    tuple val(meta), path(bam), path (bai)
    val(threshold)
    path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('pmdtools'), eval("pmdtools --version | sed 's/.*v//'")  , topic: versions, emit: versions_pmdtools
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def split_cpus = Math.floor(task.cpus/2)
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") error "[pmdtools/filter] Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    //threshold and header flags activate filtering function of pmdtools
    """
    samtools \\
        calmd \\
        ${bam} \\
        ${reference} \\
        ${args} \\
        -@ ${split_cpus} \\
    | pmdtools \\
        --threshold ${threshold} \\
        --header \\
        ${args2} \\
    | samtools \\
        view \\
        ${args3} \\
        -Sb \\
        - \\
        -@ ${split_cpus} \\
        -o ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") error "[pmdtools/filter] Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    //threshold and header flags activate filtering function of pmdtools
    """
    touch ${prefix}.bam
    """
}
