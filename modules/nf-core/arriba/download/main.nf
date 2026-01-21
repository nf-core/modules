process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/27475cdcdbcc8c0ffb6b5ca8c2e6567dbe490edb96f5df4e8f01f4f95912dcd3/data' :
        'community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654' }"

    input:
    val(genome)

    output:
    path "blacklist*${genome}*.tsv.gz"    , emit: blacklist
    path "cytobands*${genome}*.tsv"       , emit: cytobands
    path "protein_domains*${genome}*.gff3", emit: protein_domains
    path "known_fusions*${genome}*.tsv.gz", emit: known_fusions
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def arriba_version = '2.5.0'
    """
    wget https://github.com/suhrig/arriba/releases/download/v${arriba_version}/arriba_v${arriba_version}.tar.gz -O arriba_v${arriba_version}.tar.gz --no-check-certificate
    tar -xzvf arriba_v${arriba_version}.tar.gz
    rm arriba_v${arriba_version}.tar.gz
    mv arriba_v${arriba_version}/database/* .
    rm -r arriba_v${arriba_version}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba_download: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    def arriba_version = '2.5.0'
    """
    echo | gzip > blacklist_hg38_GRCh38_v${arriba_version}.tsv.gz
    touch protein_domains_hg38_GRCh38_v${arriba_version}.gff3
    touch cytobands_hg38_GRCh38_v${arriba_version}.tsv
    echo | gzip > known_fusions_hg38_GRCh38_v${arriba_version}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba_download: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
