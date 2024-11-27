process COMPUTE_GENE_SCORES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e6920e60d80922852a1b19630ebe16754cf5320d:75e2c0a29159bae8a964e43ae16a45c282fdf651-0' :
        'biocontainers/mulled-v2-e6920e60d80922852a1b19630ebe16754cf5320d:75e2c0a29159bae8a964e43ae16a45c282fdf651-0' }"

    input:
    tuple val(meta), path(normalized_counts)
    path gene_score_yaml

    output:
    tuple val(meta), path("*.txt"), emit: scores_for_mqc
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    compute_gene_scores.R $gene_score_yaml $normalized_counts $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        r-singscore: \$(Rscript -e "library(singscore); cat(as.character(packageVersion('singscore')))")
        r-GSVA: \$(Rscript -e "library(GSVA); cat(as.character(packageVersion('GSVA')))")
        r-yaml: \$(Rscript -e "library(yaml); cat(as.character(packageVersion('yaml')))")
        r-FactoMineR: \$(Rscript -e "library(FactoMineR); cat(as.character(packageVersion('FactoMineR')))")
        r-stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
    END_VERSIONS
    """

    stub:
    """
    touch scores_for_mqc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        r-singscore: \$(Rscript -e "library(singscore); cat(as.character(packageVersion('singscore')))")
        r-GSVA: \$(Rscript -e "library(GSVA); cat(as.character(packageVersion('GSVA')))")
        r-yaml: \$(Rscript -e "library(yaml); cat(as.character(packageVersion('yaml')))")
        r-FactoMineR: \$(Rscript -e "library(FactoMineR); cat(as.character(packageVersion('FactoMineR')))")
        r-stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
    END_VERSIONS
    """

}
