process PLASFLOW {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.5 bioconda::plasflow=1.1.0" // conda-forge::tensorflow=1.5.0
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasflow:1.1.0--py35_0':
        'biocontainers/plasflow:1.1.0--py35_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                  , emit: tsv
    tuple val(meta), path("*_chromosomes.fasta.gz") , emit: chromosomes
    tuple val(meta), path("*_plasmids.fasta.gz")    , emit: plasmids
    tuple val(meta), path("*_unclassified.fasta.gz"), emit: unclassified
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [[ "$assembly" == *.gz ]]; then
        gunzip -c $assembly > ${prefix}.fasta
        PlasFlow.py \\
            $args \\
            --input ${prefix}.fasta \\
            --output ${prefix}.tsv
    else
        PlasFlow.py \\
            $args \\
            --input $assembly \\
            --output ${prefix}.tsv
    fi

    if [ -f ${prefix}.tsv_chromosomes.fasta ]; then
        mv ${prefix}.tsv_chromosomes.fasta ${prefix}_chromosomes.fasta
        gzip -n ${prefix}_chromosomes.fasta
    fi

    if [ -f ${prefix}.tsv_plasmids.fasta ]; then
        mv ${prefix}.tsv_plasmids.fasta ${prefix}_plasmids.fasta
        gzip -n ${prefix}_plasmids.fasta
    fi

    if [ -f ${prefix}.tsv_unclassified.fasta ]; then
        mv ${prefix}.tsv_unclassified.fasta ${prefix}_unclassified.fasta
        gzip -n ${prefix}_unclassified.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PlasFlow: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}_chromosomes.fasta.gz
    touch ${prefix}_plasmids.fasta.gz
    touch ${prefix}_unclassified.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PlasFlow: $VERSION
    END_VERSIONS
    """
}
