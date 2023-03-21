process PICARD_LIFTOVERVCF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(input_vcf)
    path dict
    path chain
    path fasta

    output:
    tuple val(meta), path("*.lifted.vcf.gz")  , emit: vcf_lifted
    tuple val(meta), path("*.unlifted.vcf.gz"), emit: vcf_unlifted
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard LiftoverVcf] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        LiftoverVcf \\
        $args \\
        --INPUT $input_vcf \\
        --OUTPUT ${prefix}.lifted.vcf.gz \\
        --CHAIN $chain \\
        --REJECT ${prefix}.unlifted.vcf.gz \\
        --REFERENCE_SEQUENCE $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard LiftoverVcf --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lifted.vcf.gz
    touch ${prefix}.unlifted.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard LiftoverVcf --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
