process MMSEQS_CREATETAXDB {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

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
