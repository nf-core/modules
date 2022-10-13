process GATK4_COMPOSESTRTABLEFILE {
    tag "$fasta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        "${params.docker_registry ?: 'quay.io/biocontainers'}/gatk4:4.2.6.1--hdfd78af_0" }

    input:
    path(fasta)
    path(fasta_fai)
    path(dict)

    output:
    path "*.zip"            , emit: str_table
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def avail_mem = 6
    if (!task.memory) {
        log.info '[GATK ComposeSTRTableFile] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" ComposeSTRTableFile \\
        --reference $fasta \\
        --output ${fasta.baseName}.zip \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch test.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
