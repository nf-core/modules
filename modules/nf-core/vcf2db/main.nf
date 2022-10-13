process VCF2DB {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::vcf2db=2020.02.24" : null)
        'https://depot.galaxyproject.org/singularity/vcf2db:2020.02.24--hdfd78af_1':

    input:
    tuple val(meta), path(vcf), path(ped)

    output:
    tuple val(meta), path("*.db") , emit: db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2020.02.24" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    vcf2db.py \\
        $vcf \\
        $ped \\
        ${prefix}.db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2db: $VERSION
    END_VERSIONS
    """
}
