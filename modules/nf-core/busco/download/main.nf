process BUSCO_DOWNLOAD {
    tag "$lineage"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.8.2--pyhdfd78af_0':
        'biocontainers/busco:5.8.2--pyhdfd78af_0' }"

    input:
    val lineage

    output:
    path "busco_downloads", emit: download_dir
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${lineage}"
    """
    busco \\
        --download $lineage \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
