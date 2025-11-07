process EMU_ABUNDANCE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89f51e5b1961b27efc33dbafe25ae7f85c1ccfc2e2df5341849237a6db2023a1/data'
        : 'community.wave.seqera.io/library/emu:3.5.4--c9ac9572d0d77ae9'}"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}_rel-abundance.tsv")                , emit: report
    tuple val(meta), path("${prefix}_read-assignment-distributions.tsv"), emit: assignment_report, optional: true
    tuple val(meta), path("${prefix}_emu_alignments.sam")               , emit: samfile          , optional: true
    tuple val(meta), path("${prefix}_unclassified_mapped.*")            , emit: unclassified     , optional: true
    tuple val(meta), path("${prefix}_unmapped.*")                       , emit: unmapped         , optional: true
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    emu \\
        abundance \\
        ${args} \\
        --threads ${task.cpus} \\
        --db ${db} \\
        --output-basename ${prefix} \\
        ${reads}

    if [ -d results ]; then
        mv results/* .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}_rel-abundance.tsv
    touch ${prefix}_read-assignment-distributions.tsv
    touch ${prefix}_emu_alignments.sam
    touch ${prefix}_unclassified_mapped.fasta
    touch ${prefix}_unmapped.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """
}
