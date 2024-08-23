process PINTS_CALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // NOTE Stopped publishing at 1.1.9 https://quay.io/repository/biocontainers/pypints?tab=tags
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pints_caller:650ef63c9768bb5a' :
        'community.wave.seqera.io/library/pints_caller:b1e6404513917871' }"

    input:
    tuple val(meta), path(bams), path(bais)
    val assay_type

    output:
    tuple val(meta), path("*_divergent_peaks.bed")     , optional:true, emit: divergent_TREs
    tuple val(meta), path("*_bidirectional_peaks.bed") , optional:true, emit: bidirectional_TREs
    tuple val(meta), path("*_unidirectional_peaks.bed"), optional:true, emit: unidirectional_TREs
    tuple val(meta), path("peakcalling_*.log")                        , emit: peakcalling_log
    path  "versions.yml"                                              , emit: versions

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
        --bam-file $bams \\
        --save-to . \\
        --file-prefix $prefix \\
        --thread $task.cpus \\
        --dont-check-updates \\
        --exp-type $assay_type \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pints: \$(pints_caller --version)
    END_VERSIONS
    """
}
