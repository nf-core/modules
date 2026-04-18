process SOURMASH_SKETCH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sourmash:4.9.4--hdfd78af_0'
        : 'biocontainers/sourmash:4.9.4--hdfd78af_0'}"

    input:
    tuple val(meta), path(library, stageAs: 'library/*')

    output:
    tuple val(meta), path("*.sig"), emit: signatures
    tuple val("${task.process}"), val('sourmash'), eval("sourmash --version 2>&1 | sed 's/^sourmash //'"), emit: versions_sourmash, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // required defaults for the tool to run, but can be overridden
    def args = task.ext.args ?: "dna --param-string 'scaled=1000,k=21,k=31,k=51,abund'"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find -L library/ -type f > library.txt

    sourmash sketch \\
        ${args} \\
        --merge '${prefix}' \\
        --output '${prefix}.sig' \\
        --from-file library.txt
    """

    stub:
    def args = task.ext.args ?: "dna --param-string 'scaled=1000,k=21,k=31,k=51,abund'"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find -L library/ -type f > library.txt

    echo "sourmash sketch \\
        ${args} \\
        --merge '${prefix}' \\
        --output '${prefix}.sig' \\
        --from-file library.txt"

    touch '${prefix}.sig'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
