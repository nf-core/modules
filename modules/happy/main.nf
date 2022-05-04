def VERSION_HAP = '0.3.14'
def VERSION_PRE = '0.3.14'

process HAPPY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hap.py=0.3.14" : null)
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
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pre.py \\
        $args \\
        -R $bed \\
        --reference $fasta \\
        --threads $task.cpus \\
        $query_vcf \\
        ${prefix}_preprocessed.vcf.gz

    hap.py \\
        $truth_vcf \\
        ${prefix}_preprocessed.vcf.gz \\
        $args2 \\
        --reference $fasta \\
        --threads $task.cpus \\
        -R $bed \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION_HAP
        pre.py: $VERSION_PRE
    END_VERSIONS
    """
}
