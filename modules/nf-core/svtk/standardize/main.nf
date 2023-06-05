process SVTK_STANDARDIZE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::svtk=0.0.20190615"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtk:0.0.20190615--py37h73a75cf_2':
        'biocontainers/svtk:0.0.20190615--py37h73a75cf_2' }"

    input:
    tuple val(meta), path(vcf)
    path fasta_fai

    output:
    tuple val(meta), path("*.std.vcf.gz"), emit: standardized_vcf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def arguments   = args.args     ?: ''
    def caller      = args.caller   ?: 'delly'

    def VERSION = '0.0.20190615' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def contigs = fasta_fai ? "--contigs ${fasta_fai}" : ""

    """
    svtk standardize \\
        ${arguments} \\
        ${contigs} \\
        ${vcf} \\
        ${prefix}.std.vcf.gz \\
        ${caller}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtk: ${VERSION}
    END_VERSIONS
    """
}
