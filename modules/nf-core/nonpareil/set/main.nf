process NONPAREIL_SET {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::nonpareil=3.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.4.1--r42h9f5acd7_2':
        'biocontainers/nonpareil:3.4.1--r42h9f5acd7_2' }"

    input:
    tuple val(meta), path(npos)

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
    Nonpareil.set(list.files(pattern='*.npo')${args_cmd})
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
    "${task.process}": Rscript -e 'library('Nonpareil'); cat(paste(unlist(packageVersion("Nonpareil")),collapse="."))'
        ""
    END_VERSIONS
    """
}
