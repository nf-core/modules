process GATK4_LEFTALIGNANDTRIMVARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(intervals)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "versions.yml"		         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals       ? "--intervals $intervals" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK LeftAlignAndTrimVariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    gatk --java-options "-Xmx${avail_mem}G" LeftAlignAndTrimVariants \\
        $interval_command \\
        --variant $vcf \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
