process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val("${task.process}"), val('picard'), eval("picard MarkDuplicates --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${reads.getExtension()}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    if ("${reads}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    # Detect UMI (RX tag) presence in first 1000 reads
    has_umi=\$(picard ViewSam --INPUT "${reads}" | grep -v '^@' | head -n 1000 | awk 'BEGIN{n=0} /\\tRX:Z:/{n++} END{print (n>0)?1:0}')
    barcode_tag=""
    if [ "\${has_umi}" = "1" ]; then
        barcode_tag="--BARCODE_TAG RX"
    fi

    # Create index only if output will be coordinate-sorted
    # (skip if args request queryname sort order, which cannot be indexed)
    create_index=""
    if ! echo "${args}" | grep -qiE "QUERYNAME"; then
        create_index="--CREATE_INDEX true"
    fi

    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        ${args} \\
        --INPUT ${reads} \\
        --OUTPUT ${prefix}.${suffix} \\
        ${reference} \\
        \${create_index} \\
        \${barcode_tag} \\
        --METRICS_FILE ${prefix}.metrics.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${reads.getExtension()}"
    if ("${reads}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.${suffix}.bai
    touch ${prefix}.metrics.txt
    """
}
