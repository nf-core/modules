process IPHOP_DOWNLOAD {
    label 'process_single'

    conda "bioconda::iphop=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iphop:1.3.1--pyhdfd78af_0':
        'biocontainers/iphop:1.3.1--pyhdfd78af_0' }"

    input:
    val db_version

    output:
    path "iphop_db/"        , emit: iphop_db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p download_dir

    iphop \\
        download \\
        --db_dir download_dir \\
        --db_version $db_version \\
        --no_prompt \\
        $args

    iphop \\
        download \\
        --db_dir download_dir/*/ \\
        --db_version $db_version \\
        --no_prompt \\
        --full_verify

    mkdir -p iphop_db
    rm download_dir/*.tar.*
    mv download_dir/* iphop_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(iphop --version 2>&1) | head -n 1 |sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """

    """

}
