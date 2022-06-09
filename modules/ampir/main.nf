process AMPIR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-ampir=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ampir:1.1.0':
        'quay.io/biocontainers/r-ampir:1.1.0' }"

    input:
    tuple val(meta), path(faa)
    val cut_off
    val model

    output:
    tuple val(meta), path("*.faa"), emit: amps_faa
    tuple val(meta), path("*.csv"), emit: amps_csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(ampir)

    protein_seqs <- read_faa('${faa}')
    prediction <- predict_amps(protein_seqs, model = '${model}')
    prediction <- protein_seqs[which(prediction\$prob_AMP >= as.integer(${cut_off})), ]
    df_to_faa(protein_seqs, "${prefix}.faa")
    write.table(prediction, file = "${prefix}.csv", row.names = FALSE, quote = FALSE, dec = '.')

    version_file_path <- "versions.yml"
    version_ampir <- paste(unlist(packageVersion("ampir")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    ampir: ", f, sep = "")
    writeLines(version_ampir, f)
    close(f)
    """
}
