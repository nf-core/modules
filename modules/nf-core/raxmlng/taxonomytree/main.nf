process RAXMLNG_TAXONOMYTREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // quay.io/biocontainers/python:3.11 has no build-hash-suffixed tag to pin to (unlike
    // real bioconda-recipe images); pin by digest instead so the underlying image can't
    // silently drift between runs.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python@sha256:b322907f8e52b2055ccad4e46848d28a4a5631b403116cc80ddf61ec8601e05e' }"

    input:
    tuple val(meta), path(taxonomy)

    output:
    tuple val(meta), path("*.guide.nwk"), emit: guide_tree
    path "versions.yml"                 , emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'taxonomytree.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.guide.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
