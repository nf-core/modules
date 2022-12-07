process MINIPROT_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::miniprot=0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miniprot:0.5--h7132678_0':
        'quay.io/biocontainers/miniprot:0.5--h7132678_0' }"

    input:
    tuple val(meta), path(pep)
    tuple val(meta2), path(ref)

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
        miniprot: \$(miniprot --version 2>&1)
    END_VERSIONS
    """
}
