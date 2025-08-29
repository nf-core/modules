process PHARMCAT_VCFPREPROCESSOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pgkb/pharmcat:3.0.1' :
        'docker.io/pgkb/pharmcat:3.0.1' }"

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("*.vcf.bgz"), emit: preprocessor_bgz
    tuple val(meta), path("*.missing_pgx_var.vcf"), emit: missing_pgx_var
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Run preprocessor
    pharmcat_vcf_preprocessor \\
        -vcf ${vcf_file} \\
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
    # Create minimal VCF header file
    echo -e "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO" > ${prefix}.vcf

    # Duplicate .vcf file for --missing-pgx-var output
    cp ${prefix}.vcf ${prefix}.missing_pgx_var.vcf

    # Create compressed VCF stub with fallback if bgzip unavailable
    if command -v bgzip >/dev/null 2>&1; then
        bgzip -c ${prefix}.vcf > ${prefix}.vcf.bgz
    else
        gzip -c ${prefix}.vcf > ${prefix}.vcf.bgz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat --version)
    END_VERSIONS
    """
}

