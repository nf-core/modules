process POOLSNP {
    tag "$meta.id"
    label 'process_medium'
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poolsnp:1.0.0--py312h7e72e81_0':
        'biocontainers/poolsnp:1.0.0--py312h7e72e81_0' }"

    input:
    tuple val(meta)      , path(mpileup)
    tuple val(meta2)     , path(reference)
    val(max_cov)
    val(min_cov)
    val(min_count)
    val(min_freq)
    val(miss_frac)
    val(allsites)
    val(badsites)

    output:
    tuple val(meta)      , path("*.vcf.gz")   , emit: vcf
    tuple val(meta)      , path("*cov-*.txt") , emit: max_cov
    tuple val(meta)      , path("*BS.txt.gz") , emit: bad_sites, optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    PoolSNP.sh \\
        mpileup=\$PWD/${mpileup} \\
        reference=\$PWD/${reference} \\
        names=${prefix} \\
        jobs=${task.cpus} \\
        max-cov=${max_cov} \\
        min-cov=${min_cov} \\
        min-count=${min_count} \\
        min-freq=${min_freq} \\
        miss-frac=${miss_frac} \\
        badsites=${badsites} \\
        allsites=${allsites} \\
        output=\$PWD/${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        version: "${VERSION}"
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz
    touch ${prefix}_cov-${max_cov}.txt
    touch ${prefix}_BS.txt.gz

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        version: "${VERSION}"
    END_VERSIONS
    """
}
