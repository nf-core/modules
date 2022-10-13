process GATK4_CALIBRATEDRAGSTRMODEL {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container { workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        "${params.docker_registry ?: 'quay.io/biocontainers'}/gatk4:4.2.6.1--hdfd78af_0" }

    input:
    tuple val(meta), path(bam), path(bam_index), path(intervals)
    path  fasta
    path  fasta_fai
    path  dict
    path  strtablefile

    output:
    tuple val(meta), path("*.txt")   , emit: dragstr_model
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK CalibrateDragstrModel] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" CalibrateDragstrModel \\
        --input ${bam} \\
        --output ${prefix}.txt \\
        --reference ${fasta} \\
        --str-table-path ${strtablefile} \\
        --threads ${task.cpus} \\
        ${intervals_command} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
