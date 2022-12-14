process HAPPY_HAPPY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::hap.py=0.3.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.14--py27h5c5a3ab_0':
        'quay.io/biocontainers/hap.py:0.3.14--py27h5c5a3ab_0' }"

    input:
    tuple val(meta), path(truth_vcf), path(query_vcf), path(bed)
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path('*.csv'), path('*.json')  , emit: metrics
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    hap.py \\
        $truth_vcf \\
        $query_vcf \\
        $args \\
        --reference $fasta \\
        --threads $task.cpus \\
        -R $bed \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """
}
