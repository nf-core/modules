process NONPAREIL_CURVE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.5.5--r43hdcf5f25_0':
        'biocontainers/nonpareil:3.5.5--r43hdcf5f25_0' }"

    input:
    tuple val(meta), path(npo)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_cmd = args != '' ? ", ${args}" : ""
    """
    #!/usr/bin/env Rscript
    library(Nonpareil)

    png(file='${prefix}.png')
    Nonpareil.curve('${npo}'${args_cmd})
    dev.off()

    version_file_path <- "versions.yml"
    version_nonpareil <- paste(unlist(packageVersion("Nonpareil")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    nonpareil: ", f, sep = "")
    writeLines(version_nonpareil, f)
    close(f)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}": \$(Rscript -e 'library('Nonpareil'); cat(paste(unlist(packageVersion("Nonpareil")),collapse="."))')
    END_VERSIONS
    """
}
