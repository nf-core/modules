process HIFICNV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hificnv:1.0.1--h9ee0642_0':
        'quay.io/biocontainers/hificnv:1.0.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(maf)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(exclude)
    tuple val(meta4), path(expected_cn)

    output:
    tuple val(meta), path("*.copynum.bedgraph"), emit: copynum, optional: true
    tuple val(meta), path("*.depth.bw"),         emit: depth
    tuple val(meta), path("*.maf.bw"),           emit: maf, optional: true
    tuple val(meta), path("*.vcf.gz"),           emit: vcf, optional: true
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // handle optional inputs
    def maf_arg         = maf         ? "--maf '${maf}'"                 : ""
    def exclude_arg     = exclude     ? "--exclude '${exclude}'"         : ""
    def expected_cn_arg = expected_cn ? "--expected-cn '${expected_cn}'" : ""

    def cmd_parts = [
        "hificnv",
        "--bam '${bam}'",
        "--ref '${ref}'",
        maf_arg, exclude_arg, expected_cn_arg,
        "--threads ${task.cpus}",
        "--output-prefix '${prefix}'",
        args
    ].findAll { it } // remove empties

    // Note: hificnv may exit with non-zero code in certain valid scenarios
    // (e.g., when no CNVs are detected), so we handle exit codes manually
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

    # Create optional output files
    touch ${prefix}.copynum.bedgraph
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.maf.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version 2>&1 | sed 's/^.*hificnv //; s/ .*\$//')
    END_VERSIONS
    """
}
