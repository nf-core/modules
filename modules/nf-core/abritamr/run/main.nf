process ABRITAMR_RUN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6078836553c48db86c8a1d126ca764bc812b11aaeb7222299fe7be3a06ed68e/data'
        : 'community.wave.seqera.io/library/abritamr:6aee07e42fbc9d88' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.summary_matches.txt")  , emit: matches
    tuple val(meta), path("*.summary_partials.txt") , emit: partials
    tuple val(meta), path("*.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("*.amrfinder.out")        , emit: out
    tuple val(meta), path("*.abritamr.txt")         , emit: txt, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name    = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    abritamr run \\
        --contigs ${fasta_name} \\
        --prefix ${prefix} \\
        --jobs ${task.cpus} \\
        ${args}

    # Rename output files to prevent name collisions
    mv ${prefix}/summary_matches.txt ${prefix}.summary_matches.txt
    mv ${prefix}/summary_partials.txt ${prefix}.summary_partials.txt
    mv ${prefix}/summary_virulence.txt ${prefix}.summary_virulence.txt
    mv ${prefix}/amrfinder.out ${prefix}.amrfinder.out
    if [ -f results/abritamr.txt ]; then
        # This file is not always present
        mv ${prefix}/abritamr.txt ${prefix}.abritamr.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr \\([0-9.]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary_matches.txt
    touch ${prefix}.summary_partials.txt
    touch ${prefix}.summary_virulence.txt
    touch ${prefix}.amrfinder.out
    touch ${prefix}.amrfinder.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr \\([0-9.]*\\).*/\\1/')
    END_VERSIONS
    """
}
