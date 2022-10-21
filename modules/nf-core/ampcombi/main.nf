process AMPCOMBI {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ampcombi=0.1.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.1.4--pyhdfd78af_0':
        'quay.io/biocontainers/ampcombi:0.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(amp_input)
    path(faa_folder)

    output:
    tuple val(meta), path("*/*")                       , emit: results_dir
    tuple val(meta), path("*/*diamond_matches.txt")    , emit: txt
    tuple val(meta), path("*/*ampcombi.csv")           , emit: csv
    tuple val(meta), path("*/*amp.faa")                , emit: faa
    path("AMPcombi_summary.csv")                                , optional:true, emit: summary
    path("*.log")                                               , optional:true, emit: log
    path("*/amp_ref_database")                                  , optional:true, emit: results_db
    path("*/amp_ref_database/*.dmnd")                           , optional:true, emit: results_db_dmnd
    path("*/amp_ref_database/*.clean.fasta")                    , optional:true, emit: results_db_fasta
    path("*/amp_ref_database/*.tsv")                            , optional:true, emit: results_db_tsv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
//    def fileoption = amp_input.isDirectory() ? "--amp_results $amp_input" : "--path_list ${amp_input} --sample_list ${prefix}"
//    def fileoption = amp_input.isCollectionOrArray([]) ? "--path_list ${amp_input} --sample_list ${prefix}" : "--amp_results $amp_input"
//    def fileoption = amp_input.contains(['']) ? "--path_list $amp_input --sample_list ${prefix}" : "--amp_results $amp_input"
//    def fileoption = amp_input instanceof List ? "--path_list '${amp_input} --sample_list ${prefix}" : "--amp_results $amp_input"
//    def fileoption = amp_input instanceof List ? "--path_list '${amp_input.collect { "'$it'" }.split(',')}' --sample_list ${prefix}" : "--amp_results $amp_input"
//    def fileoption = amp_input instanceof List ? "--path_list ${amp_input.collect{"$it"}.toString().split(',')} --sample_list ${prefix}" : "--amp_results $amp_input"
    def fileoption = amp_input instanceof List ? "--path_list '${amp_input.collect{"$it"}.join("' '")}' --sample_list ${prefix}" : "--amp_results $amp_input"
//"[${file_list.collect { "'$it'" }.join(', ')}]"
    """
    ampcombi \\
        $args \\
        ${fileoption} \\
        --faa_folder $faa_folder/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
