process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.4.0--h0033a41_2' :
        'biocontainers/arriba:2.4.0--h0033a41_2' }"

    input:

    output:
    path "*"              , emit: reference
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz -O arriba_v2.4.0.tar.gz
    tar -xzvf arriba_v2.4.0.tar.gz
    rm arriba_v2.4.0.tar.gz
    mv arriba_v2.4.0/database/* .
    rm -r arriba_v2.4.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba_download: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_hg38_GRCh38_v2.4.0.tsv.gz
    touch protein_domains_hg38_GRCh38_v2.4.0.gff3
    touch cytobands_hg38_GRCh38_v2.4.0.tsv
    touch known_fusions_hg38_GRCh38_v2.4.0.tsv.gz
    touch protein_domains_hg38_GRCh38_v2.4.0.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba_download: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
