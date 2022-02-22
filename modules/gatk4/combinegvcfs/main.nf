process GATK4_COMBINEGVCFS {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(fastaIndex), path (fastaDict)
    path(vcffiles)
    path(vcf_idx)

    output:
    path("*.combined.g.vcf.gz"), emit: combined_gvcf
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem       = 3
    if (!task.memory) {
        log.info '[GATK COMBINEGVCFS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def variant_string = ''
    vcffiles.each {variant_string += "-V ${it} "} // loop to create a string adding -V to each vcf file
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        CombineGVCFs \\
        -R ${fasta} \\
        -O ${prefix}.combined.g.vcf.gz \\
        ${args} \\
        --tmp-dir . \\
        ${variant_string}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
