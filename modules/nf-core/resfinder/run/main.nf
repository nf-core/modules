process RESFINDER_RUN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.1.11--hdfd78af_0':
        'biocontainers/resfinder:4.1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq), path(fasta)
    path db_point
    path db_res

    output:
    tuple val(meta), path("*.json")                           , emit: json
    tuple val(meta), path("disinfinder_kma")                  , optional: true, emit: disinfinder_kma
    tuple val(meta), path("pheno_table_species.txt")          , optional: true, emit: pheno_table_species
    tuple val(meta), path("pheno_table.txt")                  , optional: true, emit: pheno_table
    tuple val(meta), path("pointfinder_kma")                  , optional: true, emit: pointfinder_kma
    tuple val(meta), path("PointFinder_prediction.txt")       , optional: true, emit: pointfinder_prediction
    tuple val(meta), path("PointFinder_results.txt")          , optional: true, emit: pointfinder_results
    tuple val(meta), path("PointFinder_table.txt")            , optional: true, emit: pointfinder_table
    tuple val(meta), path("ResFinder_Hit_in_genome_seq.fsa")  , optional: true, emit: resfinder_hit_in_genome_seq
    tuple val(meta), path("resfinder_blast")                  , optional: true, emit: resfinder_blast
    tuple val(meta), path("resfinder_kma")                    , optional: true, emit: resfinder_kma
    tuple val(meta), path("ResFinder_Resistance_gene_seq.fsa"), optional: true, emit: resfinder_resistance_gene_seq
    tuple val(meta), path("ResFinder_results_table.txt")      , optional: true, emit: resfinder_results_table
    tuple val(meta), path("ResFinder_results_tab.txt")        , optional: true, emit: resfinder_results_tab
    tuple val(meta), path("ResFinder_results.txt")            , optional: true, emit: resfinder_results
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input = ""
    if (fastq) {
        input = "-ifq " + fastq.join(" ")
    } else if (fasta) {
        input = "-ifa ${fasta}"
    }

    def db = ""
    if (db_res) {
        db = "-db_res ${db_res}"
    }
    if (db_point) {
        db = "$db -db_point ${db_point}"
    }
    """
    run_resfinder.py \\
        $args \\
        $input \\
        $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$(run_resfinder.py -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir resfinder_kma
    touch ${prefix}.json \\
        pheno_table.txt \\
        ResFinder_Hit_in_genome_seq.fsa \\
        ResFinder_Resistance_gene_seq.fsa \\
        ResFinder_results_table.txt \\
        ResFinder_results_tab.txt \\
        ResFinder_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$(run_resfinder.py -v)
    END_VERSIONS
    """
}
