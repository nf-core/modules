
process GATK4_COMBINEGVCFS {
    //tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    path(fasta)
    path(vcf)

    output:
    path("combined.gvcf.gz")      , emit: combined_gvcf
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    def avail_mem       = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        CombineGVCFs \\
        -R $fasta \\
        -O combined.gvcf.gz \\
        $args \\
        --tmp-dir . \\
        $vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
