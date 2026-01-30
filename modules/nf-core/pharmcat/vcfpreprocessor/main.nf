process PHARMCAT_VCFPREPROCESSOR {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pgkb/pharmcat:3.1.1' :
        'docker.io/pgkb/pharmcat:3.1.1' }"

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

