process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data' :
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"
    input:
    val(genome)

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
    wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz -O arriba_v2.4.0.tar.gz --no-check-certificate
    tar -xzvf arriba_v2.4.0.tar.gz
    rm arriba_v2.4.0.tar.gz
    mv arriba_v2.4.0/database/* .
    rm -r arriba_v2.4.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | sed 's/GNU Wget //' | sed -E 's/(.*) built on darwin.*/\1/')
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_hg38_GRCh38_v2.4.0.tsv.gz
    touch cytobands_hg38_GRCh38_v2.4.0.tsv
    touch protein_domains_hg38_GRCh38_v2.4.0.gff3
    touch known_fusions_hg38_GRCh38_v2.4.0.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | sed 's/GNU Wget //' | sed -E 's/(.*) built on darwin.*/\1/')
    END_VERSIONS
    """
}
