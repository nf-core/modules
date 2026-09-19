process BIGSLICE_DOWNLOADDB {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0'
        : 'quay.io/biocontainers/bigslice:2.0.2--pyh8ed023e_0'}"

    input:
    val meta

    output:
    tuple val(meta), path ("bigslice-models")              , emit: db
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('bigslice'), val("2.0.2"), topic: versions, emit: versions_bigslice
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Copy the script to the work directory so it writes bigslice-models/ here
    # (download_bigslice_hmmdb uses __file__ to determine the output path;
    # copying it here ensures bigslice-models/ is created in the work directory)

    cp \$(which download_bigslice_hmmdb) ./download_bigslice_hmmdb_local

    python \\
        ./download_bigslice_hmmdb_local \\
        ${args}

    """

    stub:
    """
    mkdir -p bigslice-models/biosyn_pfam bigslice-models/sub_pfams
    touch bigslice-models/biosyn_pfam/Biosyn_pfams.hmm
    touch bigslice-models/sub_pfams/corepfam.tsv
    """
}
