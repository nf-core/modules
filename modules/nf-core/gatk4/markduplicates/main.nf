process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*cram"),     emit: cram,  optional: true
    tuple val(meta), path("*bam"),      emit: bam,   optional: true
    tuple val(meta), path("*.crai"),    emit: crai,  optional: true
    tuple val(meta), path("*.bai"),     emit: bai,   optional: true
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = bam.collect{"--INPUT $it"}.join(' ')
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \\
        $input_list \\
        --OUTPUT ${prefix} \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        ${reference} \\
        $args


    if  [[ ${prefix} == *.cram ]]; then
        mv ${prefix}.bai ${prefix}.crai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
