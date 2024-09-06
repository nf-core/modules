process GUBBINS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.3.5--py39pl5321he4a0461_0' :
        'biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0' }"

    input:
    path alignment

    output:
    path "*.fasta"                          , emit: fasta
    path "*.gff"                            , emit: gff
    path "*.vcf"                            , emit: vcf
    path "*.csv"                            , emit: stats
    path "*.phylip"                         , emit: phylip
    path "*.recombination_predictions.embl" , emit: embl_predicted
    path "*.branch_base_reconstruction.embl", emit: embl_branch
    path "*.final_tree.tre"                 , emit: tree
    path "*.node_labelled.final_tree.tre"   , emit: tree_labelled
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir numba_cache_dir
    export NUMBA_CACHE_DIR='./numba_cache_dir'

    run_gubbins.py \\
        --threads $task.cpus \\
        $args \\
        $alignment
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir numba_cache_dir
    export NUMBA_CACHE_DIR='./numba_cache_dir'

    touch ${alignment.baseName}.fasta
    touch ${alignment.baseName}.gff
    touch ${alignment.baseName}.vcf
    touch ${alignment.baseName}.csv
    touch ${alignment.baseName}.phylip
    touch ${alignment.baseName}.recombination_predictions.embl
    touch ${alignment.baseName}.branch_base_reconstruction.embl
    touch ${alignment.baseName}.final_tree.tre
    touch ${alignment.baseName}.node_labelled.final_tree.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """
}
