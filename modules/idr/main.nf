// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IDR {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::idr=2.0.4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39hcbe4a3b_5"
    } else {
        container "quay.io/biocontainers/idr:2.0.4.2--py38h9af456f_5"
    }

    input:
    path peaks
    val peak_type
    val prefix

    output:
    path "*idrValues.txt", emit: idr
    path "*log.txt"      , emit: log
    path "*.png"         , emit: png
    path "versions.yml"  , emit: versions

    script:
    // TODO: Check if number of peaks is 1 and generate error or fail if so
        // https://github.com/nf-core/modules/blob/f47c27edfbc6181779699ff6b919d773d578aed7/modules/picard/mergesamfiles/main.nf#L37
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[IDR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
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
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(idr --version 2>&1) | sed 's/^.*IDR //; s/ .*\$//')
    END_VERSIONS
    """
}
