process KRONA_KTUPDATETAXONOMY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a763169a99fb3b4ccea1102edab08c60fc2888f42852f7dd2540c80434c504c/data'
        : 'community.wave.seqera.io/library/krona_make:7bb1fe2561793909'}"

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    tuple val("${task.process}"), val('krona'), eval("ktImportTaxonomy | grep -Po '(?<=KronaTools )[0-9.]+'"), topic: versions, emit: versions_krona

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ktUpdateTaxonomy.sh \\
        ${args} \\
        taxonomy/
    """

    stub:
    """
    mkdir taxonomy
    touch taxonomy/taxonomy.tab
    """
}
