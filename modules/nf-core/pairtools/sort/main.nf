process PAIRTOOLS_SORT {
    tag "$meta.id"
    label 'process_high'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "bioconda::pairtools=1.0.2 conda-forge::numpy=1.23"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.0.2--py39h2a9f597_0' :
        'biocontainers/pairtools:1.0.2--py39h2a9f597_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: sorted
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem = (task.memory.giga*0.8).intValue()
    """
    pairtools \\
        sort \\
        $args \\
        --nproc $task.cpus \\
        --memory "${mem}G" \\
        -o ${prefix}.pairs.gz \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
