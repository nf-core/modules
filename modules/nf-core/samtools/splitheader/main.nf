process SAMTOOLS_SPLITHEADER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_readgroups.txt"), emit: readgroup
    tuple val(meta), path("*_programs.txt"), emit: programs
    tuple val(meta), path("*_sequences.txt"), emit: sequences
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        -H \\
        ${args} \\
        ${input} \\
    | tee \\
        >( grep '^@RG' > ${prefix}_readgroups.txt ) \
        >( grep '^@PG' > ${prefix}_programs.txt ) \
        >( grep '^@SQ' > ${prefix}_sequences.txt ) \
     > /dev/null
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA" > ${prefix}_readgroups.txt
    echo -e "@PG\\tID:samtools.4\\tPN:samtools\\tPP:samtools.3\\tVN:1.22.1\\tCL:samtools view -H ${input}" > ${prefix}_programs.txt
    echo -e "@SQ\\tSN:chr1\\tLN:10000" > ${prefix}_sequences.txt
    """
}
