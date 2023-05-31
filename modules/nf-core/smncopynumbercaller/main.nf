process SMNCOPYNUMBERCALLER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::smncopynumbercaller=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smncopynumbercaller:1.1.2--py310h7cba7a3_0' :
        'biocontainers/smncopynumbercaller:1.1.2--py310h7cba7a3_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("out/*.tsv"),  emit: smncopynumber
    tuple val(meta), path("out/*.json"), emit: run_metrics
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.1.2" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    echo $bam | tr ' ' '
    ' > manifest.txt
    smn_caller.py \\
        $args \\
        --manifest manifest.txt \\
        --prefix $prefix \\
        --outDir "out" \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SMNCopyNumberCaller: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.1.2"
    """
    mkdir out
    touch out/${prefix}.tsv
    touch out/${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SMNCopyNumberCaller: $VERSION
    END_VERSIONS
    """
}
