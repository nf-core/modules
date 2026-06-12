process NONPAREIL_CURVE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.5.5--r43hdcf5f25_0':
        'quay.io/biocontainers/nonpareil:3.5.5--r43hdcf5f25_0' }"

    input:
    tuple val(meta), path(npo)

    output:
    tuple val(meta), path("*.png"), optional: true, emit: png
    tuple val("${task.process}"), val('Nonpareil'), eval('Rscript -e "library(Nonpareil); cat(paste(unlist(packageVersion(\'Nonpareil\')), collapse = \'.\'))"'), emit: versions_nonpareil, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_cmd = args != '' ? ", ${args}" : ""
    """
    Rscript - <<EOF
    library(Nonpareil)
    png(file='${prefix}.png')
    Nonpareil.curve('${npo}'${args_cmd})
    dev.off()
    EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    """
}
