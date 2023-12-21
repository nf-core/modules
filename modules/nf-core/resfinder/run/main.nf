process RESFINDER_RUN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.1.11--hdfd78af_0':
        'biocontainers/resfinder:4.1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq), path(fasta)
    path db

    output:
    tuple val(meta), path("*.json")                           , emit: json
    tuple val(meta), path("disinfinder_kma")                  , emit: disinfinder_kma
    tuple val(meta), path("pointfinder_kma")                  , emit: pointfinder_kma
    tuple val(meta), path("pheno_table.txt")                  , emit: pheno_table
    tuple val(meta), path("ResFinder_Hit_in_genome_seq.fsa")  , emit: hit_in_genome_seq
    tuple val(meta), path("resfinder_kma")                    , emit: resfinder_kma
    tuple val(meta), path("ResFinder_Resistance_gene_seq.fsa"), emit: resistance_gene_seq
    tuple val(meta), path("ResFinder_results_table.txt")      , emit: results_table
    tuple val(meta), path("ResFinder_results_tab.txt")        , emit: results_tab
    tuple val(meta), path("ResFinder_results.txt")            , emit: results
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
    """
    run_resfinder.py \\
        $args \\
        $input \\
        -d $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$(run_resfinder.py -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir disinfinder_kma \\
        pointfinder_kma \\
        resfinder_kma
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
