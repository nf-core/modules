process PINTS_CALLER {
    tag "$meta.id"
    label 'process_medium'

    conda    (params.enable_conda ? "bioconda::pypints=1.1.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypints:1.1.6--pyh5e36f6f_1' :
        'quay.io/biocontainers/pypints:1.1.6--pyh5e36f6f_1' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*_divergent_peaks.bed")     , emit: divergent_TREs
    tuple val(meta), path("*_bidirectional_peaks.bed") , emit: bidirectional_TREs
    tuple val(meta), path("*_unidirectional_peaks.bed"), emit: unidirectional_TREs
    path  "versions.yml"     , emit: versions

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
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pints: \$(pints_caller --version)
    END_VERSIONS
    """
}
