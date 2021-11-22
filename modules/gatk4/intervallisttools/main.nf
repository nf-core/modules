process GATK4_INTERVALLISTTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(interval_list)

    output:
    tuple val(meta), path("*_split/*/*.interval_list"), emit: interval_list
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """

    mkdir ${prefix}_split

    gatk \\
        IntervalListTools \\
        -I ${interval_list} \\
        -O ${prefix}_split \\
        $args

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
    ${task.process.tokenize(':').last()}:
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
