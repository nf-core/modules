
process MINIPROT_ALIGN {

    tag "$meta.id"
    label 'process_medium'

    def version = '0.5-c2'
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the miniprot process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/miniprot:${version}"

    input:
    tuple val(meta), path(pep)
    path ref


    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--gff") ? "gff" : "paf"
    """
    miniprot \\
        $args \\
        -t $task.cpus \\
        ${ref} \\
        ${pep} \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: $version
    END_VERSIONS
    """
}

