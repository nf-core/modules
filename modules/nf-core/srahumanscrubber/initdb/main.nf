process SRAHUMANSCRUBBER_INITDB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.0.0--hdfd78af_0':
        'quay.io/biocontainers/sra-human-scrubber:2.0.0--hdfd78af_0' }"

    output:
    path "*.human_filter.db", emit: db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '2.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    DBVERSION=\$(curl "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/current/version.txt")
    curl -f "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/\${DBVERSION}.human_filter.db" -o "\${DBVERSION}.human_filter.db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: $VERSION
        sra-human-scrubber-db: \$DBVERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '2.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    DBVERSION="STUBRUN"
    touch "\${DBVERSION}.human_filter.db"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: $VERSION
        sra-human-scrubber-db: $VERSION
    END_VERSIONS
    """
}
