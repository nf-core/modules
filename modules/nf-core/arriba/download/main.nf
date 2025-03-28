process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.4.0--h0033a41_2' :
        'biocontainers/arriba:2.4.0--h0033a41_2' }"

    input:
    val(genome)
    val(arriba_database)

    output:
    path "blacklist*${genome}*.tsv.gz"       , emit: blacklist
    path "cytobands*${genome}*.tsv"          , emit: cytobands
    path "protein_domains*${genome}*.gff3"   , emit: protein_domains
    path "known_fusions*${genome}*.tsv.gz"   , emit: known_fusions
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget $arriba_database -O arriba_database.tar.gz --no-check-certificate
    mkdir arriba_database
    ## Ensures --strip-components only applied when top level of tar contents is a directory
    ## If just files or multiple directories, place all in prefix
    if [[ \$(tar -taf arriba_database.tar.gz | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -xavf \\
            arriba_database.tar.gz \\
            -C arriba_database --strip-components 1
    else
        tar \\
            -xavf \\
            arriba_database.tar.gz \\
            -C arriba_database
    fi
    rm arriba_database.tar.gz
    mv arriba_database/database/* .
    rm -r arriba_database

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba_database: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_hg38_GRCh38_v2.4.0.tsv.gz
    touch protein_domains_hg38_GRCh38_v2.4.0.gff3
    touch cytobands_hg38_GRCh38_v2.4.0.tsv
    touch known_fusions_hg38_GRCh38_v2.4.0.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba_atabase: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
