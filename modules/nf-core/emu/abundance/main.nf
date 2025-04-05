process EMU_ABUNDANCE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emu:3.5.1--hdfd78af_0':
        'biocontainers/emu:3.5.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.tsv")                                     , emit: report
    tuple val(meta), path("*read-assignment-distributions.tsv")        , emit: assignment_report, optional:true
    tuple val(meta), path("*.sam")                                     , emit: samfile, optional:true
    tuple val(meta), path("*.fasta")                                   , emit: unclassified_fa , optional:true
    path "versions.yml"                                                , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    emu \\
        abundance \\
        $args \\
        --threads $task.cpus \\
        --db $db \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.sam
    touch ${prefix}.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """
}
