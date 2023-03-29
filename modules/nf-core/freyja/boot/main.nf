process FREYJA_BOOT {
    tag "$meta.id"
    label 'process_long'

    conda "bioconda::freyja=1.3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.3.12--pyhdfd78af_0':
        'quay.io/biocontainers/freyja:1.3.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(variants), path(depths)
    val repeats
    path barcodes
    path lineages_meta

    output:
    tuple val(meta), path("*lineages.csv")  , emit: lineages
    tuple val(meta), path("*summarized.csv"), emit: summarized
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freyja \\
        boot \\
        $args \\
        --nt $task.cpus \\
        --nb $repeats \\
        --output_base $prefix \\
        --barcodes $barcodes \\
        --meta $lineages_meta \\
        $variants \\
        $depths

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_lineage.csv
    touch ${prefix}_summarized.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
