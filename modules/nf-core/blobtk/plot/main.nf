process BLOBTK_PLOT {
    // Linked to issue https://github.com/sanger-tol/genomenote/issues/184
    // Depending on the blob dataset in use, the grid option may not
    // work at all. This is down to the version of blobtoolkit used to
    // generate the blob.
    // Adding a check would overly complicate the module so for now
    // we can ignore errors, with the knowledge it would only kill
    // runs in which the blobdir doesn't have the right data.
    errorStrategy 'ignore'

    tag "$prefix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/243d043f1c9e152e75dbb0ef8c64022df50efbcaa4e1bbaea36bebd751e84e93/data' :
        'community.wave.seqera.io/library/blobtk:0.7.1--e3f63bb2cdc8fb96' }"

    input:
    tuple val(meta), path(fasta)
    path(local_path)                // Genuine path location must be a path.
    val(online_path)                // HTTPS location needs to remain a value
    val extra_args                  // In format [name: "", args: ""]

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''

    if ( online_path && local_path ) {
        error "BLOBTK_PLOT can't use both local_path and online_path, use `[]` as input for the unused channel."
    }

    def resource = online_path ?: local_path
    prefix       = task.ext.prefix ?: "${meta.id}"

    """
    blobtk plot \\
        -d $resource \\
        $args \\
        -o ${prefix}.png
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtk: \$(blobtk --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
