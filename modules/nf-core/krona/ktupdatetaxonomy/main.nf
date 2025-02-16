process KRONA_KTUPDATETAXONOMY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a763169a99fb3b4ccea1102edab08c60fc2888f42852f7dd2540c80434c504c/data':
        'community.wave.seqera.io/library/krona_make:7bb1fe2561793909' }"

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    ktUpdateTaxonomy.sh \\
        $args \\
        taxonomy/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$(ktImportTaxonomy | grep -Po "(?<=KronaTools )[0-9.]+")
    END_VERSIONS
    """

    stub:
    """
    mkdir taxonomy

    touch taxonomy/taxonomy.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$(ktImportTaxonomy | grep -Po "(?<=KronaTools )[0-9.]+")
    END_VERSIONS
    """
}
