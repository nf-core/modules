process PRIDEPY_FETCHSDRF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pridepy:0.0.15--pyhdfd78af_0':
        'quay.io/biocontainers/pridepy:0.0.15--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(pride_id)

    output:
    tuple val(meta), path("${prefix}.sdrf.tsv"), emit: sdrf
    tuple val("${task.process}"), val('pridepy'), eval('python -c "from importlib.metadata import version; print(version(\'pridepy\'))"'), topic: versions, emit: versions_pridepy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pridepy stream-files-metadata \\
        -a "${pride_id}" \\
        -o files_metadata.json \\
        $args

    # Check if SDRF is present in metadata and extract
    sdrf_name=\$(awk -F'"' '/"fileName" *: *"(sdrf\\.tsv|[^"]+\\.sdrf\\.tsv)"/ {print \$4; exit}' files_metadata.json)
    [ -n "\$sdrf_name" ] || { echo "ERROR: No SDRF file found for ${pride_id}" >&2; exit 1; }

    pridepy download-file-by-name \\
        -a "${pride_id}" \\
        -f "\$sdrf_name" \\
        -o . \\
        $args

    mv "\$sdrf_name" "${prefix}.sdrf.tsv"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sdrf.tsv
    """
}
