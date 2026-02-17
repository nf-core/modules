process HAPPY_PREPY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.15--py27hcb73b3d_0':
        'biocontainers/hap.py:0.3.15--py27hcb73b3d_0' }"

    input:
    tuple val(meta), path(vcf), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path('*.vcf.gz')  , emit: preprocessed_vcf
    tuple val("${task.process}"), val('happy'), val('0.3.15'), topic: versions, emit: versions_happy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def restrict_region = bed ? "-R ${bed}": ""
    """
    pre.py \\
        $args \\
        $restrict_region \\
        --reference $fasta \\
        --threads $task.cpus \\
        $vcf \\
        ${prefix}.vcf.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz

    """
}
