process ENSEMBLVEP_FILTERVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3da6e21cbf9803529421d7e136d1ebec5ff71ec50e0d996eda2ce11ec2c19bf9/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf:1e13f65f931a6954'}"

    input:
    tuple val(meta), path(input)
    path feature_file

    output:
    tuple val(meta), path("*.${extension}"), emit: output
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep
    tuple val("${task.process}"), val('perl-math-cdf'), eval("perl -MMath::CDF -e 'print \\\$Math::CDF::VERSION'"), topic: versions, emit: versions_perlmathcdf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.suffix ?: "vcf"
    """
    filter_vep \\
        ${args} \\
        --input_file ${input} \\
        --output_file ${prefix}.${extension} \\
        --only_matched
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.suffix ?: "vcf"
    """
    touch ${prefix}.${extension}
    """
}
