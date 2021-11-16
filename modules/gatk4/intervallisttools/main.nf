// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_INTERVALLISTTOOLS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(interval_list)

    output:
    tuple val(meta), path("*_split/*/*.interval_list"), emit: interval_list
    path "versions.yml"           , emit: versions

    script:
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

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
