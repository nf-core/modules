process PRIDEPY_DOWNLOADFILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pridepy:0.0.15--pyhdfd78af_0':
        'quay.io/biocontainers/pridepy:0.0.15--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(file_name), val(pride_accession)

    output:
    tuple val(meta), path("${file_name}"), emit: file
    tuple val("${task.process}"), val('pridepy'), eval('python -c "from importlib.metadata import version; print(version(\'pridepy\'))"'), topic: versions, emit: versions_pridepy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    pridepy download-file-by-name \\
        -a "${pride_accession}" \\
        -f "${file_name}" \\
        -o . \\
        ${args}

    # pridepy exits 0 even when the download fails, so fail explicitly on an empty file
    [ -s "${file_name}" ] || { echo "ERROR: Downloaded file ${file_name} is empty" >&2; exit 1; }
    """

    stub:
    """
    touch "${file_name}"
    """
}
