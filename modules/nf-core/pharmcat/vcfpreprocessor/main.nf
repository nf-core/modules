process PHARMCAT_VCFPREPROCESSOR {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a/data' :
        'community.wave.seqera.io/library/pharmcat3:3.1.1--876b7152770ba008' }"

    input:
    tuple val(meta), path(vcf_file), path(vcf_idx)

    output:
    tuple val(meta), path("*.vcf.bgz") 		   , emit: preprocessor_bgz
    tuple val(meta), path("*.missing_pgx_var.vcf") , emit: missing_pgx_var
    path "versions.yml"				   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pharmcat_vcf_preprocessor \\
        -vcf ${vcf_file} \\
        --concurrent-mode \\
        --max-concurrent-processes $task.cpus \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.missing_pgx_var.vcf
    echo | gzip > ${prefix}.vcf.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat --version)
    END_VERSIONS
    """
}
