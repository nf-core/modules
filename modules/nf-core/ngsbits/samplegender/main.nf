process NGSBITS_SAMPLEGENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngs-bits:2024_11--py312hd80e9a6_0':
        'biocontainers/ngs-bits:2024_11--py312hd80e9a6_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    val method

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref = fasta ? "-ref ${fasta}" : ""
    """
    SampleGender \\
        -in ${bam} \\
        -method ${method} \\
        -out ${prefix}.tsv \\
        ${ref} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(SampleGender --version 2>&1) | sed 's/SampleGender //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(SampleGender --version 2>&1) | sed 's/SampleGender //' )
    END_VERSIONS
    """
}
