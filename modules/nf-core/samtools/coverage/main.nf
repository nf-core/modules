process SAMTOOLS_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.txt"), emit: coverage
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    if (input.name.endsWith('.cram') && (!fasta || !fai)) {
        error("CRAM input file provided but no reference FASTA and/or FAI index for said reference, both are required for CRAM input.")
    }
    """
    samtools \\
        coverage \\
        ${args} \\
        -o ${prefix}.txt \\
        ${reference} \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > ${prefix}.txt
    """
}
