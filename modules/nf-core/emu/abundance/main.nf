process EMU_ABUNDANCE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/emu:3.5.1--hdfd78af_0'
        : 'biocontainers/emu:3.5.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}_rel-abundance.tsv")                , emit: report
    tuple val(meta), path("${prefix}_read-assignment-distributions.tsv"), emit: assignment_report, optional: true
    tuple val(meta), path("${prefix}_emu_alignments.sam")               , emit: samfile          , optional: true
    tuple val(meta), path("${prefix}_unclassified_mapped.fasta")        , emit: unclassified_fa  , optional: true
    tuple val(meta), path("${prefix}_unmapped.fasta")                   , emit: unmapped_fa      , optional: true
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
