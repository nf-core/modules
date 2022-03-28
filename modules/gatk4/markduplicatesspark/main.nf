process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'broadinstitute/gatk:4.2.3.0' }"

    input:
    tuple val(meta), path(bams)
    path  fasta
    path  fasta_fai
    path  dict

    output:
    tuple val(meta), path("${prefix}"), emit: output
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bam_list = bams.collect{"--INPUT ${it}"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    export SPARK_USER=spark3

    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicatesSpark \\
        $bam_list \\
        --reference ${fasta} \\
        --tmp-dir . \\
        --output ${prefix} \\
        --spark-master local[${task.cpus}] \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
