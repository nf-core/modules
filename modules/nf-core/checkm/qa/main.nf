process CHECKM_QA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6e77f70239b60110da040c4307b8048749ed1fc86262e07d27f1eb12a314d14f/data'
        : 'community.wave.seqera.io/library/checkm-genome:1.2.5--8d1d1a2477a013ce'}"

    input:
    tuple val(meta), path(analysis_dir), path(marker_file), path(coverage_file)
    path exclude_marker_file
    path db

    output:
    tuple val(meta), path("${prefix}.txt"), optional: true, emit: output
    tuple val(meta), path("${prefix}.fasta"), optional: true, emit: fasta
    tuple val("${task.process}"), val('checkm'), eval("checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//'"), emit: versions_checkm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def checkm_db = db ? "export CHECKM_DATA_PATH=${db}" : ""
    suffix = task.ext.args?.matches(".*-o 9.*|.*--out_file 9.*") ? "fasta" : "txt"
    def coverage = coverage_file && coverage_file.isFile() ? "--coverage_file ${coverage_file}" : ""
    def exclude = exclude_marker_file && exclude_marker_file.isFile() ? "--exclude_markers ${exclude_marker_file}" : ""
    """
    ${checkm_db}

    checkm \\
        qa \\
        --threads ${task.cpus} \\
        --file ${prefix}.${suffix} \\
        ${marker_file} \\
        ${analysis_dir} \\
        ${coverage} \\
        ${exclude} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p checkm_dummy_db
    export CHECKM_DATA_PATH=\$PWD/checkm_dummy_db
    touch ${prefix}.txt ${prefix}.fasta
    """
}
