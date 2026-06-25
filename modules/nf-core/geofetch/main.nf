process GEOFETCH {
    tag "$geo_accession"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/geofetch:0.12.6--pyh7cba7a3_0':
        'quay.io/biocontainers/geofetch:0.12.6--pyh7cba7a3_0' }"

    input:
    val geo_accession

    output:
    tuple val("${geo_accession}"), path("${geo_accession}/*.CEL.gz"), emit: samples
    tuple val("${task.process}"), val('geofetch'), eval("geofetch --version 2>&1 | sed '1!d; s/^geofetch //'"), topic: versions, emit: versions_geofetch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    geofetch \\
        -i \\
        $geo_accession \\
        --processed \\
        -g . \\
        $args
    """

    stub:

    """
    mkdir -p ${geo_accession}
    echo "" | gzip > ${geo_accession}/foo.CEL.gz
    """
}
