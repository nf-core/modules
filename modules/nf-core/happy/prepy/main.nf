process HAPPY_PREPY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::hap.py=0.3.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.14--py27h5c5a3ab_0':
        'biocontainers/hap.py:0.3.14--py27h5c5a3ab_0' }"

    input:
    tuple val(meta), path(vcf), path(bed)
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path('*.vcf.gz')  , emit: preprocessed_vcf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    pre.py \\
        $args \\
        -R $bed \\
        --reference $fasta \\
        --threads $task.cpus \\
        $vcf \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pre.py: $VERSION
    END_VERSIONS
    """
}
