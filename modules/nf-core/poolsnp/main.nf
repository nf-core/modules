process POOLSNP {
    tag "$meta.id"
    label 'process_medium'
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poolsnp:1.0.1--py312h7e72e81_0':
        'biocontainers/poolsnp:1.0.1--py312h7e72e81_0' }"

    input:
    tuple val(meta) , path(mpileup)
    tuple val(meta2), path(reference)
    tuple val(meta) , val(max_cov), path(max_cov_file)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    tuple val(meta), path("*cov-*.txt"), emit: max_cov  , optional: true
    tuple val(meta), path("*BS.txt.gz"), emit: bad_sites, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    assert (!max_cov && max_cov_file) || (max_cov && !max_cov_file)

    """
    PoolSNP.sh \\
        mpileup=\$PWD/${mpileup} \\
        output=\$PWD/${prefix} \\
        names=${prefix} \\
        reference=\$PWD/${reference} \\
        jobs=${task.cpus} \\
        max-cov=${max_cov ? "${max_cov}" : "\$PWD/${max_cov_file}"} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        poolsnp: "${VERSION}"
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo "##fileformat=VCFv4.2" > ${prefix}.vcf
    echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${prefix}.vcf
    gzip ${prefix}.vcf
    ${max_cov ? "touch ${prefix}_cov-${max_cov}.txt" : ""}
    echo "" | gzip > ${prefix}_BS.txt.gz

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        poolsnp: "${VERSION}"
    END_VERSIONS
    """
}
