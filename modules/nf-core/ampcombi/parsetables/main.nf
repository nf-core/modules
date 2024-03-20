process AMPCOMBI_PARSETABLES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.2.1--pyhdfd78af_0':
        'biocontainers/ampcombi:0.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(amp_input)
    path( faa_input )
    path( gbk_input )
    path( opt_amp_db )

    output:
    tuple val(meta), path("${meta.id}*")                        , emit: sample_dir
    tuple val(meta), path("${meta.id}/*diamond_matches.txt")    , emit: txt
    tuple val(meta), path("${meta.id}/*ampcombi.tsv")           , emit: tsv
    tuple val(meta), path("${meta.id}/*amp.faa")                , emit: faa
    tuple val(meta), path("${meta.id}/*_ampcombi.log")          , optional:true, emit: sample_log
    tuple val(meta), path("Ampcombi_parse_tables.log")          , optional:true, emit: full_log
    tuple val(meta), path("*/amp_ref_database")                 , optional:true, emit: results_db
    tuple val(meta), path("*/amp_ref_database/*.dmnd")          , optional:true, emit: results_db_dmnd
    tuple val(meta), path("*/amp_ref_database/*.clean.fasta")   , optional:true, emit: results_db_fasta
    tuple val(meta), path("*/amp_ref_database/*.tsv")           , optional:true, emit: results_db_tsv
    path "versions.yml"                                         , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db = opt_amp_db? "--amp_database $opt_amp_db": ""
    def faa = faa_input.isDirectory() ? "--faa ${faa_input}/" : "--faa ${faa_input}"
    def gbk = gbk_input.isDirectory() ? "--gbk ${gbk_input}/" : "--gbk ${gbk_input}"
    """
    ampcombi parse_tables \\
        $args \\
        --path_list '${amp_input.collect{"$it"}.join("' '")}' \\
        --sample_list ${prefix} \\
        --aminoacid_length 2000 \\
        --db_evalue 2000 \\
        --amp_cutoff 0 \\
        --ampir_file '.tsv' \\
        --amplify_file '.tsv' \\
        --macrel_file '.prediction' \\
        --neubi_file '.fasta' \\
        --hmmsearch_file 'candidates.txt' \\
        --ampgram_file '.tsv' \\
        --amptransformer_file '.txt' \\
        --threads $task.cpus \\
        ${db} \\
        ${faa} \\
        ${gbk}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi parse_tables --version | sed 's/ampcombi//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db = opt_amp_db? "--amp_database $opt_amp_db": ""
    """
    mkdir -p ${prefix}
    touch ${prefix.id}/*diamond_matches.txt
    touch ${prefix}/*ampcombi.tsv
    touch ${prefix}/*amp.faa
    touch ${prefix}/*_ampcombi.log
    touch Ampcombi_parse_tables.log

    mkdir -p amp_ref_database
    touch amp_ref_database/*.dmnd
    touch amp_ref_database/*.clean.fasta
    touch amp_ref_database/*.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi//')
    END_VERSIONS
    """
}
