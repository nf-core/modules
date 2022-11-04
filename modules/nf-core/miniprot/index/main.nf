process MINIPROT_INDEX {
    tag "$meta.id"
    label 'process_medium'
    def version = '0.5-c2'
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the miniprot process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/miniprot:${version}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mpi"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    miniprot \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mpi \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: $version
    END_VERSIONS
    """
}
