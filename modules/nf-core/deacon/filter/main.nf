// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEACON_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f1270fe75323ac6d213bf66d72b9f292c5ff31f7f5ee12300d2bed2b858f5ab/data'
        : 'community.wave.seqera.io/library/deacon:0.20.0--09d7e8e5c3f9dc88'}"

    input:
    tuple val(meta), path(index), path(reads)

    output:
    tuple val(meta), path(fastq_out)                , emit: fastq_filtered
    tuple val(meta), path("${prefix}.json")         , emit: log
    tuple val(meta), path("${prefix}.debug.tsv.gz") , emit: debug, optional: true
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    fastq_out = reads instanceof List ? ["${prefix}_1.fq.gz", "${prefix}_2.fq.gz"] : ["${prefix}.fq.gz"]
    if ([reads].flatten().any { read -> "${read}" in fastq_out*.toString() }) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    def output_reads = reads instanceof List ? "-o ${fastq_out[0]} -O ${fastq_out[1]}" : "-o ${fastq_out[0]}"
    // --debug needs a path: the module names it so that the file is declared as an output
    def debug_flag = args.contains('--debug') ? "--debug ${prefix}.debug.tsv.gz" : ''
    def args_corrected = args.replaceAll(/--debug(?:[ =]+[^-\s]\S*)?/, '').trim()
    """
    deacon \\
        filter \\
        --threads ${task.cpus} \\
        ${args_corrected} \\
        ${debug_flag} \\
        --summary ${prefix}.json \\
        -d $index \\
        $reads \\
        ${output_reads}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    fastq_out = reads instanceof List ? ["${prefix}_1.fq.gz", "${prefix}_2.fq.gz"] : ["${prefix}.fq.gz"]
    def debug_stub = args.contains('--debug') ? "echo | gzip > ${prefix}.debug.tsv.gz" : ''
    """
    for fastq in ${fastq_out.join(' ')}; do
        echo | gzip > \$fastq
    done
    touch ${prefix}.json
    ${debug_stub}
    """
}
