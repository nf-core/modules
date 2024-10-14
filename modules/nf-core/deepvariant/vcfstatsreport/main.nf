process DEEPVARIANT_VCFSTATSREPORT {
    tag "$meta.id"
    label 'process_single'

    // FIXME Conda is not supported at the moment
    // BUG https://github.com/nf-core/modules/issues/1754
    // BUG https://github.com/bioconda/bioconda-recipes/issues/30310
    container "nf-core/deepvariant:1.6.1"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.visual_report.html"), emit: report
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // WARN https://github.com/nf-core/modules/pull/5801#issuecomment-2194293755
    // FIXME Revert this on next version bump
    def VERSION = '1.6.1'

    """
    /opt/deepvariant/bin/vcf_stats_report \\
        --input_vcf=${vcf} \\
        --outfile_base ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    // WARN https://github.com/nf-core/modules/pull/5801#issuecomment-2194293755
    // FIXME Revert this on next version bump
    def VERSION = '1.6.1'
    """
    touch ${prefix}.visual_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: $VERSION
    END_VERSIONS
    """
}
