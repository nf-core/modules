process SPOLIGOTYPER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/18627f630e2f90d14aa45effad5909fb77092cb48534b7d677f91f029c65d2de/data'
        : 'community.wave.seqera.io/library/spoligotyper:0.9.0--679a58efe8efff70'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_spoligotyping.txt")    , emit: tsv
    tuple val(meta), path("${prefix}_spoligotyping.json")   , emit: json
    tuple val(meta), path("${prefix}_spoligotyping.pdf")    , emit: pdf, optional: true
    tuple val(meta), path("${prefix}_spoligotyping_mqc.json"), emit: multiqc
    tuple val("${task.process}"), val('spoligotyper'), eval("spoligotyper --version | sed 's/spoligotyper //'"), emit: versions_spoligotyper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Reads (single-end or paired-end) or an assembly: spoligotyper detects the file type
    def input = reads instanceof List && reads.size() == 2 ? "-r1 \"${reads[0]}\" -r2 \"${reads[1]}\"" : "-r1 \"${reads}\""
    """
    spoligotyper \\
        $args \\
        $input \\
        --sample "${prefix}" \\
        --output . \\
        --threads $task.cpus
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_spoligotyping.txt" "${prefix}_spoligotyping.json" "${prefix}_spoligotyping.pdf" "${prefix}_spoligotyping_mqc.json"
    """
}
