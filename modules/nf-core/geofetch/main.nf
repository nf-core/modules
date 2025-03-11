process GEOFETCH {
    tag "$geo_accession"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/geofetch:0.12.6--pyh7cba7a3_0':
        'biocontainers/geofetch:0.12.6--pyh7cba7a3_0' }"

    input:
    val geo_accession

    output:
    tuple val("${geo_accession}"), path("${geo_accession}/*.CEL.gz"), emit: samples
    path "versions.yml"                                             , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        geofetch: \$(geofetch --version|& sed '1!d ; s/geofetch //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir -p ${geo_accession}
    cd ${geo_accession}
    touch foo.CEL
    gzip foo.CEL
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        geofetch: \$(geofetch --version|& sed '1!d ; s/geofetch //')
    END_VERSIONS
    """
}
