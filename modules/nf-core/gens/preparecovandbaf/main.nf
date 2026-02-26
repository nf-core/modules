process PREPARECOVANDBAF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/67abc0d3d1eaaeaff3eafb36babecf56b5912c2e7b0c5100a9f32eda1c82cb30/data':
        'community.wave.seqera.io/library/htslib_python_pip_gens-input-data-tools:8fd1a0ecd4a60110' }"

    input:
    tuple val(meta), path(read_counts), path(gvcf), path(gvcf_tbi)
    path baf_positions

    output:
    tuple val(meta), path("*.cov.bed.gz")     , emit: cov_gz
    tuple val(meta), path("*.cov.bed.gz.tbi") , emit: cov_tbi
    tuple val(meta), path("*.baf.bed.gz")     , emit: baf_gz
    tuple val(meta), path("*.baf.bed.gz.tbi") , emit: baf_tbi
    tuple val("${task.process}"), val('preparecovandbaf'), eval("generate_cov_and_baf --version"), topic: versions, emit: versions_preparecovandbaf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_cov_and_baf \\
        --coverage $read_counts \\
        --gvcf $gvcf \\
        --label $prefix \\
        --baf_positions $baf_positions \\
        --bgzip_tabix_output \\
        $args \\
        --outdir .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.cov.bed.gz
    touch ${prefix}.cov.bed.gz.tbi
    echo "" | gzip > ${prefix}.baf.bed.gz
    touch ${prefix}.baf.bed.gz.tbi
    """
}
