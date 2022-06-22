VERSION = "2020.02.24"

process VCF2DB {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcf2db=2020.02.24" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf2db:2020.02.24--hdfd78af_1':
        'quay.io/biocontainers/vcf2db:2020.02.24--hdfd78af_1' }"

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
