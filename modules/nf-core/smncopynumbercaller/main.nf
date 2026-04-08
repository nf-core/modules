process SMNCOPYNUMBERCALLER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smncopynumbercaller:1.1.2--py310h7cba7a3_0' :
        'biocontainers/smncopynumbercaller:1.1.2--py310h7cba7a3_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("out/*.tsv"),  emit: smncopynumber
    tuple val(meta), path("out/*.json"), emit: run_metrics
    tuple val("${task.process}"), val('smncopynumbercaller'), val('1.1.2'), topic: versions, emit: versions_smncopynumbercaller
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $bam | tr ' ' '
    ' > manifest.txt
    smn_caller.py \\
        $args \\
        --manifest manifest.txt \\
        --prefix $prefix \\
        --outDir "out" \\
        --threads $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir out
    touch out/${prefix}.tsv
    touch out/${prefix}.json
    """
}
