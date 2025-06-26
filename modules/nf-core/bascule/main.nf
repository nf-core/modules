process BASCULE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c46a35446e5fc41efe84d4962b1900fdb7dda26333fb0a44ecb7c521dcda1ba6/data':
        'community.wave.seqera.io/library/pybascule_r-bascule_python_r-ggplot2_pruned:797049b34cb704fc' }"

    input:
    tuple val(meta), path(counts_matrices)  // counts_matrices = folder with .csv files

    output:
    tuple val(meta), path("*_fit_bascule.rds")  , emit: bascule_rds
    tuple val(meta), path("*_plots_bascule.rds"), emit: bascule_plots_rds
    tuple val(meta), path("*_plots_bascule.pdf"), emit: bascule_plots_pdf
    tuple val(meta), path("*_plots_bascule.png"), emit: bascule_plots_png
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}_fit_bascule.rds
    touch ${prefix}_plots_bascule.rds
    touch ${prefix}_plots_bascule.pdf
    touch ${prefix}_plots_bascule.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bascule: \$(Rscript -e 'library(bascule); sessionInfo()\$otherPkgs\$bascule\$Version')
        ggplot2: \$(Rscript -e 'library(ggplot2); sessionInfo()\$otherPkgs\$ggplot2\$Version')
        reticulate: \$(Rscript -e 'library(reticulate); sessionInfo()\$otherPkgs\$reticulate\$Version')
        tidyverse: \$(Rscript -e 'library(tidyverse); sessionInfo()\$otherPkgs\$tidyverse\$Version')
    END_VERSIONS
    """
}
