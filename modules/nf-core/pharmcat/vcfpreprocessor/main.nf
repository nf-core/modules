process PHARMCAT_VCFPREPROCESSOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.1.1--876b7152770ba008'}"

    input:
    tuple val(meta), path(vcf_gz), path(vcf_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(pharmcat_positions), path(pharmcat_positions_index)
    tuple val(meta5), path(pharmcat_uniallelic_positions), path(pharmcat_uniallelic_positions_index)

    output:
    tuple val(meta), path("*.preprocessed.vcf.bgz"),                                                                                                    emit: preprocessed_vcf
    tuple val(meta), path("*.missing_pgx_var.vcf"),                                                                                 optional: true,     emit: missing_pgx_var
    tuple val("${task.process}"), val('pharmcat_vcf_preprocessor'), eval("pharmcat_vcf_preprocessor --version | cut -f4 -d ' '"),   topic: versions,    emit: versions_pharmcat_vcf_preprocessor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pharmcat_vcf_preprocessor \\
        --vcf ${vcf_gz} \\
        --base-filename ${prefix} \\
        --reference-genome ${fasta} \\
        --reference-pgx-vcf ${pharmcat_positions} \\
        --output-dir . \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.preprocessed.vcf.bgz
    touch ${prefix}.missing_pgx_var.vcf
    """
}
