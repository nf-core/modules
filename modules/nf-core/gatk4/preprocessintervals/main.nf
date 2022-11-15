process GATK4_PREPROCESSINTERVALS {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    path  exclude_intervals
    path  fasta
    path  fai
    path  dict

    output:
    path  "*.interval_list"                 , emit: interval_list
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exclude_command = exclude_intervals ? "--exclude-intervals $exclude_intervals" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK PreprocessIntervals] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" PreprocessIntervals \\
        $exclude_command \\
        --reference $fasta \\
        --padding 0 \\
        -imr OVERLAPPING_ONLY \\
        --output preprocessed_intervals.interval_list \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
