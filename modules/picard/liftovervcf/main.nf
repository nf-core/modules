process PICARD_LIFTOVERVCF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.27.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(input_vcf)
    path dict
    path chain
    path fasta

    output:
    tuple val(meta), path("*lifted.vcf.gz")  , emit: vcf_lifted
    tuple val(meta), path("*unlifted.vcf.gz"), emit: vcf_unlifted
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 1
    if (!task.memory) {
        log.info '[Picard LiftoverVcf] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
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
