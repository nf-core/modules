process GATK4_SAMTOFASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path  "versions.yml"               , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def output   = meta.single_end ? "FASTQ=${prefix}.fastq.gz" : "FASTQ=${prefix}_1.fastq.gz SECOND_END_FASTQ=${prefix}_2.fastq.gz"
    """
    gatk SamToFastq \\
        I=$bam \\
        $output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
