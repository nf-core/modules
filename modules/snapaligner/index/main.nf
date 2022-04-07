process SNAPALIGNER_INDEX {
    tag '$fasta'
    label 'process_high'

    conda (params.enable_conda ? "bioconda::snap-aligner=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap-aligner:2.0.1--hd03093a_1':
        'quay.io/biocontainers/snap-aligner:2.0.1--hd03093a_1' }"

    input:
    path fasta

    output:
    path "snap"            ,emit: index
    path "versions.yml"    ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir snap
    snap-aligner \\
    index \\
    $fasta \\
    snap \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //')
    END_VERSIONS
    """
    stub:
    """
    mkdir snap
    echo "Genome" > snap/Genome
    echo "GenomeIndex" > snap/GenomeIndex
    echo "GenomeIndexHash" > snap/GenomeIndexHash
    echo "OverflowTable" > snap/OverflowTable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //')
    END_VERSIONS
    """
}
