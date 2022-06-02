process GATK_INDELREALIGNER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk=3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk:3.5--hdfd78af_11':
        'quay.io/biocontainers/gatk:3.5--hdfd78af_11' }"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    path path(fasta)
    path(fai)
    path(dict)
    path(known_vcf)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known = known_vcf ? "-known ${known_vcf}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK IndelRealigner] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk3 \\
        -T IndelRealigner \\
        -R ${fasta} \\
        -nt ${task.cpus}
        -I ${bam} \\
        -targetIntervals ${intervals} \\
        ${known} \\
        -o ${prefix}.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """
}
