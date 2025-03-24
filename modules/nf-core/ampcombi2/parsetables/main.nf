process AMPCOMBI2_PARSETABLES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(amp_input)
    path faa_input
    path gbk_input
    val opt_amp_db
    path opt_amp_db_dir
    path opt_interproscan

    output:
    tuple val(meta), path("${meta.id}/")                                , emit: sample_dir
    tuple val(meta), path("${meta.id}/contig_gbks/")                    , emit: contig_gbks  , optional:true
    tuple val(meta), path("${meta.id}/${meta.id}_mmseqs_matches.tsv")   , emit: db_tsv       , optional:true
    tuple val(meta), path("${meta.id}/${meta.id}_ampcombi.tsv")         , emit: tsv          , optional:true
    tuple val(meta), path("${meta.id}/${meta.id}_amp.faa")              , emit: faa          , optional:true
    tuple val(meta), path("${meta.id}/${meta.id}_ampcombi.log")         , emit: sample_log   , optional:true
    tuple val(meta), path("Ampcombi_parse_tables.log")                  , emit: full_log     , optional:true
    tuple val(meta), path("amp_${opt_amp_db}_database/")                , emit: db           , optional:true
    tuple val(meta), path("amp_${opt_amp_db}_database/*.txt")           , emit: db_txt       , optional:true
    tuple val(meta), path("amp_${opt_amp_db}_database/*.fasta")         , emit: db_fasta     , optional:true
    tuple val(meta), path("amp_${opt_amp_db}_database/mmseqs2/")        , emit: db_mmseqs    , optional:true
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_dir = opt_amp_db_dir ? "--amp_database_dir ${opt_amp_db_dir}" : ""
    def interpro = opt_interproscan ? "--interproscan_output ${opt_interproscan}" : ""

    """
    ampcombi parse_tables \\
        --path_list '${amp_input.collect { "${it}" }.join("' '")}' \\
        --faa ${faa_input} \\
        --gbk ${gbk_input} \\
        --sample_list ${prefix} \\
        --amp_database ${opt_amp_db} \\
        ${db_dir} \\
        ${interpro} \\
        ${args} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_dir = opt_amp_db_dir ? "--amp_database_dir ${opt_amp_db_dir}" : ""
    def interpro = opt_interproscan ? "--interproscan_output ${opt_interproscan}" : ""

    """
    mkdir -p ${prefix}
    mkdir -p ${prefix}/contig_gbks
    touch ${prefix}/${meta.id}_mmseqs_matches.tsv
    touch ${prefix}/${meta.id}_ampcombi.tsv
    touch ${prefix}/${meta.id}_amp.faa
    touch ${prefix}/${meta.id}_ampcombi.log
    touch Ampcombi_parse_tables.log

    mkdir -p amp_${opt_amp_db}_database
    mkdir -p amp_${opt_amp_db}_database/mmseqs2
    touch amp_${opt_amp_db}_database/*.fasta
    touch amp_${opt_amp_db}_database/*.txt
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB.dbtype
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB_h
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB_h.dbtype
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB_h.index
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB.index
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB.lookup
    touch amp_${opt_amp_db}_database/mmseqs2/ref_DB.source

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
