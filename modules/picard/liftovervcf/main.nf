process PICARD_LIFTOVERVCF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(input_vcf)
    path dict
    path chain
    path fasta

    output:
    tuple val(meta), path("*lifted.vcf")  , emit: lifted
    tuple val(meta), path("*unlifted.vcf"), emit: unlifted
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 1

    """
    picard \\
        -Xmx${avail_mem}g \\
        LiftoverVcf \\
        $args \\
        I=$input_vcf \\
        O=${prefix}.lifted.vcf \\
        CHAIN=$chain \\
        REJECT=${prefix}.unlifted.vcf \\
        R=$fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard LiftoverVcf --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
