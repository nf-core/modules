process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'quay.io/biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    tuple val(meta), path(reads)
    path  fasta
    path  index

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $reads \\
        -o ${prefix}.sam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphmap2: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
