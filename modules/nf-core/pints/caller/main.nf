process PINTS_CALLER {
    tag "$meta.id" + "${chr_name ? ' | ' + chr_name : ''}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // NOTE Stopped publishing at 1.1.9 https://quay.io/repository/biocontainers/pypints?tab=tags
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1a9e30012e1b41baf9acd1ff94e01161138d8aa17f4e97aa32f2dc4effafcd1/data' :
        'community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:39699b96998ec5f6' }"

    input:
    tuple val(meta), path(bam, arity: '1'), val(chr_name)
    val assay_type

    output:
    tuple val(meta), path("*_1_divergent_peaks.bed")     , optional:true, emit: divergent_TREs
    tuple val(meta), path("*_1_bidirectional_peaks.bed") , optional:true, emit: bidirectional_TREs
    tuple val(meta), path("*_1_unidirectional_peaks.bed"), optional:true, emit: unidirectional_TREs
    tuple val(meta), path("peakcalling_*.log")                        , emit: peakcalling_log
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" + (chr_name ? '_' + chr_name : '_all')
    // TODO handle bigwigs
    // def input_type  = ("${input[0]}".endsWith(".bam")) ? "--bam-file $input" :
    //                    ("$input".contains(".bw")) ? "--bw-pl ${input[0]} --bw-mn ${input[1]}" :
    //                    error "Please use bam or BigWig files"
    def chr = chr_name ? "--chromosome-start-with $chr_name" : ''
    """
    pints_caller \\
        --bam-file $bam \\
        --save-to . \\
        --file-prefix $prefix \\
        --thread $task.cpus \\
        --dont-check-updates \\
        --exp-type $assay_type \\
        $chr \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pints: \$(pints_caller --version)
    END_VERSIONS
    """
}
