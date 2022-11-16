process HLALA_TYPING {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::hla-la=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.3--hd03093a_0':
        'quay.io/biocontainers/hla-la:1.0.3--hd03093a_0' }"

    input:
    tuple val(meta), path(bam)
    path(graph)

    output:
    tuple val(meta), path("${prefix}")   , emit: folder
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    HLA-LA.pl \\
        --BAM $bam \\
        --graph $graph \\
        --sampleID $prefix \\
        --workingDir . \\
        --maxThreads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hla-la: 1.0.3
    END_VERSIONS
    """
}
