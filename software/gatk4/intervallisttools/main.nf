// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_INTERVALLISTTOOLS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(interval_list)

    output:
    tuple val(meta), path("*_split/*/*.interval_list"), emit: interval_list
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """

    mkdir ${prefix}_split

    gatk \\
    IntervalListTools \\
    -I ${interval_list} \\
    -O ${prefix}_split \\
    $options.args

    python3 <<CODE
    import glob, os
    # The following python code snippet rename the output files into different name to avoid overwriting or name conflict
    intervals = sorted(glob.glob("*_split/*/*.interval_list"))
    for i, interval in enumerate(intervals):
        (directory, filename) = os.path.split(interval)
        newName = os.path.join(directory, str(i + 1) + filename)
        os.rename(interval, newName)
    CODE

    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}
