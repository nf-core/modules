process FREYJA_UPDATE {
    label 'process_single'

    conda "bioconda::freyja=1.3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.3.12--pyhdfd78af_0':
        'quay.io/biocontainers/freyja:1.3.12--pyhdfd78af_0' }"

    output:
    path "${prefix}/usher_barcodes.csv"     , emit: barcodes
    path "${prefix}/curated_lineages.json"  , emit: lineages
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "freyja_db"

    """
    mkdir $prefix

    freyja \\
        update \\
        --outdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version\s+//' )
    END_VERSIONS
    """
}
