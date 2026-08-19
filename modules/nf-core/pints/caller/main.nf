process PINTS_CALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // NOTE Stopped publishing at 1.1.9 https://quay.io/repository/biocontainers/pypints?tab=tags
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1a9e30012e1b41baf9acd1ff94e01161138d8aa17f4e97aa32f2dc4effafcd1/data' :
        'community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:39699b96998ec5f6' }"

    input:
    tuple val(meta), path(bams), path(bais)
    val assay_type

    output:
    tuple val(meta), path("*_divergent_peaks.bed")     , emit: divergent_TREs     , optional:true
    tuple val(meta), path("*_bidirectional_peaks.bed") , emit: bidirectional_TREs , optional:true
    tuple val(meta), path("*_unidirectional_peaks.bed"), emit: unidirectional_TREs, optional:true
    tuple val(meta), path("peakcalling_*.log")         , emit: peakcalling_log
    tuple val("${task.process}"), val('pints') , eval("pints_caller --version")               , emit: versions_pints , topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO handle bigwigs
    // def input_type  = ("${input[0]}".endsWith(".bam")) ? "--bam-file $input" :
    //                    ("$input".contains(".bw")) ? "--bw-pl ${input[0]} --bw-mn ${input[1]}" :
    //                    error "Please use bam or BigWig files"
    """
    pints_caller \\
        --bam-file ${bams} \\
        --save-to . \\
        --file-prefix ${prefix} \\
        --thread ${task.cpus} \\
        --dont-check-updates \\
        --exp-type ${assay_type} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_divergent_peaks.bed
    touch ${prefix}_bidirectional_peaks.bed
    touch ${prefix}_unidirectional_peaks.bed
    touch peakcalling_2026_08_05_14_20_53_41.log
    """
}
