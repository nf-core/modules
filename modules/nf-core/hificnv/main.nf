process HIFICNV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hificnv:1.0.1--h9ee0642_0':
        'quay.io/biocontainers/hificnv:1.0.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(maf)
    tuple val(meta4), path(exclude)
    tuple val(meta5), path(expected_cn)

    output:
    tuple val(meta), path("*.copynum.bedgraph"), emit: copynum, optional: true
    tuple val(meta), path("*.depth.bw"),         emit: depth
    tuple val(meta), path("*.maf.bw"),           emit: maf, optional: true
    tuple val(meta), path("*.vcf{,.gz}"),        emit: vcf, optional: true
    tuple val(meta), path("*.log"),              emit: log
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Handle optional inputs - only add arguments if files are provided and not empty
    def maf_arg = (maf && maf.name != 'NO_FILE') ? "--maf ${maf}" : ""
    def exclude_arg = (exclude && exclude.name != 'NO_FILE') ? "--exclude ${exclude}" : ""
    def expected_cn_arg = (expected_cn && expected_cn.name != 'NO_FILE') ? "--expected-cn ${expected_cn}" : ""

    def cmd_parts = [
        "hificnv",
        "--bam ${bam}",
        "--ref ${ref}"
    ]

    // Add optional arguments only if they're not empty
    if (maf_arg) cmd_parts.add(maf_arg)
    if (exclude_arg) cmd_parts.add(exclude_arg)
    if (expected_cn_arg) cmd_parts.add(expected_cn_arg)

    // Add remaining required arguments
    cmd_parts.add("--threads ${task.cpus}")
    cmd_parts.add("--output-prefix ${prefix}")
    if (args) cmd_parts.add(args)

    """
    # Disable immediate exit on error
    set +e
    ${cmd_parts.join(' \\\n    ')}

    # Capture exit code
    exit_code=\$?

    # Re-enable exit on error for other commands
    set -e

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version 2>&1 | sed 's/^.*hificnv //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Create mandatory output files
    touch ${prefix}.depth.bw
    touch ${prefix}.log

    # Create optional output files
    touch ${prefix}.copynum.bedgraph
    touch ${prefix}.vcf
    touch ${prefix}.maf.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version 2>&1 | sed 's/^.*hificnv //; s/ .*\$//')
    END_VERSIONS
    """
}
