process GATK_INDELREALIGNER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk=3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk:3.5--hdfd78af_11':
        'quay.io/biocontainers/gatk:3.5--hdfd78af_11' }"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta), path(fasta)
    tuple val(meta), path(known_vcf)

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
    """
    gatk3 \\
        -T RealignerTargetCreator \\
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
