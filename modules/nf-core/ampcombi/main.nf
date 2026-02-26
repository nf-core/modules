process AMPCOMBI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.1.7--pyhdfd78af_0':
        'biocontainers/ampcombi:0.1.7--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(amp_input)
    path(faa_input)
    path(opt_amp_db)

    output:
    tuple val(meta), path("${meta.id}/")                    , emit: sample_dir
    tuple val(meta), path("${meta.id}/*diamond_matches.txt"), emit: txt
    tuple val(meta), path("${meta.id}/*ampcombi.csv")       , emit: csv
    tuple val(meta), path("${meta.id}/*amp.faa")            , emit: faa
    tuple val(meta), path("AMPcombi_summary.csv")           , emit: summary_csv     , optional:true
    tuple val(meta), path("AMPcombi_summary.html")          , emit: summary_html    , optional:true
    tuple val(meta), path("*.log")                          , emit: log             , optional:true
    tuple val(meta), path("amp_ref_database/")              , emit: results_db      , optional:true
    tuple val(meta), path("amp_ref_database/*.dmnd")        , emit: results_db_dmnd , optional:true
    tuple val(meta), path("amp_ref_database/*.clean.fasta") , emit: results_db_fasta, optional:true
    tuple val(meta), path("amp_ref_database/*.tsv")         , emit: results_db_tsv  , optional:true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated.

    Reason:
    This module is no longer recommended for use to parse results from antimicrobial tools.
    It is recommended to use ampcombi v.0.2.2 submodules instead:
    - nf-core/modules/ampcombi2/parse_tables
    - nf-core/modules/ampcombi2/complete
    - nf-core/modules/ampcombi2/cluster

    """
    assert false: deprecation_message
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db = opt_amp_db? "--amp_database $opt_amp_db": ""
    """
    ampcombi \\
        --path_list '${amp_input.collect{file_path -> "$file_path"}.join("' '")}' \\
        --sample_list ${prefix} \\
        ${db} \\
        --faa ${faa_input} \\
        ${args} \\
        --log True \\
        --threads ${task.cpus} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated.

    Reason:
    This module is no longer recommended for use to parse results from antimicrobial tools.
    It is recommended to use ampcombi v.0.2.2 submodules instead:
    - nf-core/modules/ampcombi2/parse_tables
    - nf-core/modules/ampcombi2/complete
    - nf-core/modules/ampcombi2/cluster

    """
    assert false: deprecation_message
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/*diamond_matches.txt
    touch ${prefix}/*ampcombi.csv
    touch ${prefix}/*amp.faa
    touch AMPcombi_summary.csv
    touch AMPcombi_summary.html
    touch *.log

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
