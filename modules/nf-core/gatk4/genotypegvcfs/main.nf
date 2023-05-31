process GATK4_GENOTYPEGVCFS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(gvcf), path(gvcf_index), path(intervals), path(intervals_index)
    path  fasta
    path  fai
    path  dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gvcf_command = gvcf.name.endsWith(".vcf") || gvcf.name.endsWith(".vcf.gz") ? "$gvcf" : "gendb://$gvcf"
    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    def interval_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GenotypeGVCFs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" GenotypeGVCFs \\
        --variant $gvcf_command \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        $interval_command \\
        $dbsnp_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
