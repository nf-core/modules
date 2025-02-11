process PHAROKKA_INSTALLDATABASES {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharokka:1.7.3--pyhdfd78af_0':
        'biocontainers/pharokka:1.7.3--pyhdfd78af_0' }"

    input:

    output:
    path("${prefix}/")      , emit: pharokka_db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'pharokka_db'
    """
    install_databases.py \\
        --outdir $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'pharokka_db'
    """
    mkdir -p $prefix
    touch $prefix/1Aug2023_data.tsv
    touch $prefix/1Aug2023_genomes.fa.msh
    touch $prefix/CARD
    touch $prefix/CARD.dbtype
    touch $prefix/CARD.index
    touch $prefix/CARD.lookup
    touch $prefix/CARD.source
    touch $prefix/CARD_h
    touch $prefix/CARD_h.dbtype
    touch $prefix/CARD_h.index
    touch $prefix/VFDB_setB_pro.fas.gz
    touch $prefix/VFDBclusterRes_cluster.tsv
    touch $prefix/VFDBclusterRes_rep_seq.fasta
    touch $prefix/all_phrogs.h3m
    touch $prefix/aro_index.tsv
    touch $prefix/phrog_annot_v4.tsv
    touch $prefix/phrogs_db
    touch $prefix/phrogs_db.dbtype
    touch $prefix/phrogs_db.index
    touch $prefix/phrogs_profile_db
    touch $prefix/phrogs_profile_db.dbtype
    touch $prefix/phrogs_profile_db.index
    touch $prefix/phrogs_profile_db_consensus
    touch $prefix/phrogs_profile_db_consensus.dbtype
    touch $prefix/phrogs_profile_db_consensus.index:
    touch $prefix/phrogs_profile_db_h
    touch $prefix/phrogs_profile_db_h.index
    touch $prefix/phrogs_profile_db_seq
    touch $prefix/phrogs_profile_db_seq.dbtype
    touch $prefix/phrogs_profile_db_seq.index
    touch $prefix/phrogs_profile_db_seq_h
    touch $prefix/phrogs_profile_db_seq_h.index
    touch $prefix/vfdb
    touch $prefix/vfdb.dbtype
    touch $prefix/vfdb.index
    touch $prefix/vfdb.lookup
    touch $prefix/vfdb.source
    touch $prefix/vfdb_h
    touch $prefix/vfdb_h.dbtype
    touch $prefix/vfdb_h.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """
}
