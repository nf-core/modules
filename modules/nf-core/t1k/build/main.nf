process T1K_BUILD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c64676623c43855d3a9dd8a5b02a03b13216e454a8bd1e9a7d9d5695daf5030b/data'
        : 'community.wave.seqera.io/library/t1k:1.0.9--793c68c25d680b97'}"

    input:
    tuple val(meta), path(ena), path(fasta), path(annotation)

    output:
    tuple val(meta), path("*_rna_seq.fa")  , optional: true, emit: rna_sequences
    tuple val(meta), path("*_dna_seq.fa")  , optional: true, emit: dna_sequences
    tuple val(meta), path("*_rna_coord.fa"), optional: true, emit: rna_coordinates
    tuple val(meta), path("*_dna_coord.fa"), optional: true, emit: dna_coordinates
    tuple val(meta), path("*.dat")         , optional: true, emit: database
    tuple val("${task.process}"), val('t1k'), eval("run-t1k 2>&1 | head -n 1 | cut -d '-' -f 1 | cut -d v -f 2"), topic: versions, emit: versions_t1k

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def ena_args = ena ? "-d ${ena}" : ''
    def fasta_args = fasta ? "-f ${fasta}" : ''
    def annotation_args = annotation ? "-g ${annotation}" : ''
    """
    t1k-build.pl \\
        ${args} \\
        --prefix ${prefix} \\
        ${ena_args} \\
        ${fasta_args} \\
        ${annotation_args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_rna_seq.fa
    touch ${prefix}_dna_seq.fa
    touch ${prefix}_rna_coord.fa
    touch ${prefix}_dna_coord.fa
    touch ${prefix}.dat
    """
}
