// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MMSEQS_CREATETAXDB {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/161e2ded0fc1b8d724333808695d7da12a8643668a59b652796e7a51ff9b4f3e/data'
        : 'community.wave.seqera.io/library/mmseqs2:20.ae60e--0746edfc81192f32'}"

    input:
    tuple val(meta), path(db)
    tuple val(meta2), path(taxdump_dir)
    tuple val(meta3), path(tax_mapping_file)

    output:
    tuple val(meta), path(db), emit: db_with_taxonomy
    tuple val("${task.process}"), val('mmseqs'), eval("mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //'"), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    taxdump_opt = taxdump_dir ? "--ncbi-tax-dump ${taxdump_dir}" : ""
    tax_mapping_opt = taxdump_dir && tax_mapping_file ? "--tax-mapping-file ${tax_mapping_file}" : ""
    """
    # Find database files, remove extension, and sort to get base name first (e.g., 'mmseqs' before 'mmseqs_aln')
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' | sort | head -1)

    mmseqs createtaxdb \\
      \${DB_INPUT_PATH_NAME} \\
      ./tmp \\
      --threads ${task.cpus} \\
      ${taxdump_opt} \\
      ${tax_mapping_opt} \\
      ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    """
    # Find database files, remove extension, and sort to get base name first (e.g., 'mmseqs' before 'mmseqs_aln')
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' | sort | head -1)

    echo ${args}
    touch "${db}/\$(basename \${DB_INPUT_PATH_NAME})_mapping"
    touch "${db}/\$(basename \${DB_INPUT_PATH_NAME})_taxonomy"
    """
}
