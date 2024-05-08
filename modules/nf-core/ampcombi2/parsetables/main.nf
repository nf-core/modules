process AMPCOMBI2_PARSETABLES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.2.2--pyhdfd78af_0':
        'biocontainers/ampcombi:0.2.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(amp_input)
    path(faa_input)
    path(gbk_input)
    path(opt_amp_db)

    output:
    tuple val(meta), path("${meta.id}/")                              , emit: sample_dir
    tuple val(meta), path("${meta.id}/contig_gbks/")                  , emit: contig_gbks
    tuple val(meta), path("${meta.id}/${meta.id}_diamond_matches.txt"), emit: txt
    tuple val(meta), path("${meta.id}/${meta.id}_ampcombi.tsv")       , emit: tsv
    tuple val(meta), path("${meta.id}/${meta.id}_amp.faa")            , emit: faa
    tuple val(meta), path("${meta.id}/${meta.id}_ampcombi.log")       , emit: sample_log, optional:true
    tuple val(meta), path("Ampcombi_parse_tables.log")                , emit: full_log, optional:true
    tuple val(meta), path("amp_ref_database/")                        , emit: results_db, optional:true
    tuple val(meta), path("amp_ref_database/*.dmnd")                  , emit: results_db_dmnd, optional:true
    tuple val(meta), path("amp_ref_database/*.clean.fasta")           , emit: results_db_fasta, optional:true
    tuple val(meta), path("amp_ref_database/*.tsv")                   , emit: results_db_tsv, optional:true
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db = opt_amp_db? "--amp_database $opt_amp_db": ""
    """
    ampcombi parse_tables \\
    --path_list '${amp_input.collect{"$it"}.join("' '")}' \\
    --faa ${faa_input} \\
    --gbk ${gbk_input} \\
    --sample_list ${prefix} \\
    ${db} \\
    $args \\
    --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db = opt_amp_db? "--amp_database $opt_amp_db": ""
    """
    mkdir -p ${prefix}
    mkdir -p ${prefix}/contig_gbks
    touch ${prefix}/${meta.id}_diamond_matches.txt
    touch ${prefix}/${meta.id}_ampcombi.tsv
    touch ${prefix}/${meta.id}_amp.faa
    touch ${prefix}/${meta.id}_ampcombi.log
    touch Ampcombi_parse_tables.log

    mkdir -p amp_ref_database
    touch amp_ref_database/*.dmnd
    touch amp_ref_database/*.clean.fasta
    touch amp_ref_database/*.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
