process SVTK_VCF2BED {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::svtk=0.0.20190615"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtk:0.0.20190615--py37h73a75cf_2':
        'biocontainers/svtk:0.0.20190615--py37h73a75cf_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '0.0.20190615' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    svtk vcf2bed \\
        ${vcf} \\
        ${prefix}.bed \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtk: ${VERSION}
    END_VERSIONS
    """
}
