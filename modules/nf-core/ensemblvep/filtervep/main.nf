// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ENSEMBLVEP_FILTERVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6febf703e6df28e6259f8b382a75e16cd7039be1da722a8ee716004479b82458/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf_htslib_gzip_tar:43f77a48cc79bf38'}"

    input:
    tuple val(meta), path(input)
    path feature_file
    val extension

    output:
    tuple val(meta), path("*.${extension}"), emit: output
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep
    tuple val("${task.process}"), val('perl-math-cdf'), eval("perl -MMath::CDF -e 'print \\\$Math::CDF::VERSION'"), topic: versions, emit: versions_perlmathcdf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = extension ?: "vcf"
    """
    filter_vep \\
        ${args} \\
        --input_file ${input} \\
        --output_file ${prefix}.${suffix} \\
        --only_matched
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = extension ?: "vcf"
    """
    touch ${prefix}.${suffix}
    """
}
