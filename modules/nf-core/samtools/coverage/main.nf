process SAMTOOLS_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

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
    echo "chr21\t16570000\t16610000\t8741\t39996\t99.9875\t32.4854\t29.6\t59.8" >> ${prefix}.txt
    """
}
