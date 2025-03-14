process PINTS_CALLER {
    tag "$meta.id" + "${chr_name ? ' | ' + chr_name : ''}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // NOTE Stopped publishing at 1.1.9 https://quay.io/repository/biocontainers/pypints?tab=tags
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e16d0837ec2c28760a4567b621ec4516bd780dfa8df5601eba73358a814e8c8/data' :
        'community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:47d773fa19345f15' }"

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
        pints: \$(pints_caller --version)
    END_VERSIONS
    """
}
