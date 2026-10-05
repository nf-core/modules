process SPOLIGOTYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spoligotyper:0.9.0--pyhdfd78af_0' :
        'quay.io/biocontainers/spoligotyper:0.9.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_spoligotyping.txt")    , emit: tsv
    tuple val(meta), path("*_spoligotyping.json")   , emit: json
    tuple val(meta), path("*_spoligotyping.pdf")    , emit: pdf, optional: true
    tuple val(meta), path("*_spoligotyping_mqc.json"), emit: multiqc
    tuple val("${task.process}"), val('spoligotyper'), eval("spoligotyper --version | sed 's/spoligotyper //'"), emit: versions_spoligotyper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Reads (single-end or paired-end) or an assembly: spoligotyper detects the file type
    def input = reads instanceof List && reads.size() == 2 ? "-r1 ${reads[0]} -r2 ${reads[1]}" : "-r1 ${reads}"
    """
    spoligotyper \\
        $args \\
        $input \\
        --sample ${prefix} \\
        --output . \\
        --threads $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_spoligotyping.txt ${prefix}_spoligotyping.json ${prefix}_spoligotyping.pdf ${prefix}_spoligotyping_mqc.json
    """
}
