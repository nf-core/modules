process GATK_REALIGNERTARGETCREATOR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk=3.5" : null)
    def container_image = "/gatk:3.5--hdfd78af_11"
                                                              container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(input), path(index)
    path fasta
    path fai
    path dict
    path known_vcf

    output:
    tuple val(meta), path("*.intervals"), emit: intervals
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known = known_vcf ? "-known ${known_vcf}" : ""
    if ("$input" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK RealignerTargetCreator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk3 \\
        -Xmx${avail_mem}g \\
        -T RealignerTargetCreator \\
        -nt ${task.cpus} \\
        -I ${input} \\
        -R ${fasta} \\
        -o ${prefix}.intervals \\
        ${known} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """
}
