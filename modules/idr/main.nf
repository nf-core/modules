process IDR {
    tag "$prefix"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::idr=2.0.4.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39hcbe4a3b_5' :
        'quay.io/biocontainers/idr:2.0.4.2--py39hcbe4a3b_5' }"

    input:
    path peaks
    val peak_type
    val prefix

    output:
    path "*idrValues.txt", emit: idr
    path "*log.txt"      , emit: log
    path "*.png"         , emit: png
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (peaks.toList().size < 2) {
        log.error "[ERROR] idr needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    def idr_vals = prefix ? "${prefix}.idrValues.txt" : "idrValues.txt"
    def log_file = prefix ? "${prefix}.log.txt" : "log.txt"
    """
    idr \\
        --samples $peaks \\
        --input-file-type $peak_type \\
        --output-file $idr_vals \\
        --log-output-file $log_file \\
        --plot \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        idr: \$(echo \$(idr --version 2>&1) | sed 's/^.*IDR //; s/ .*\$//')
    END_VERSIONS
    """
}
