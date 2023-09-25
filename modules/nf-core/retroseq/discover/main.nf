process RETROSEQ_DISCOVER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::perl-retroseq=1.5=pl5321hdfd78af_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-retroseq:1.5--pl5321hdfd78af_1':
        'biocontainers/perl-retroseq:1.5--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.5" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    retroseq.pl \\
        -discover \\
        -bam $bam \\
        -refTEs $tab\\
        -output ${prefix}.tab\\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retroseq_discover: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.5" 
    """
    touch out/${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retroseq_discover: $VERSION
    END_VERSIONS
    """
}