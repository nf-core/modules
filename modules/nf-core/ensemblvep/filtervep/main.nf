process ENSEMBLVEP_FILTERVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/9021c17a89ce9034b23f85736fd7ce906c44716e612f8e34d2a5e7ceeeb372d8/data'
        : 'community.wave.seqera.io/library/ensembl-vep_htslib_perl-math-cdf_gzip_tar:35e34ac4f7e9d58a'}"

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
