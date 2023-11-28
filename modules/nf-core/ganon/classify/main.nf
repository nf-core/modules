process GANON_CLASSIFY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ganon:1.5.1--py310h8abeb55_0':
        'biocontainers/ganon:1.5.1--py310h8abeb55_0' }"

    input:
    tuple val(meta) , path(fastqs)
    path(db)

    output:
    tuple val(meta), path("*.tre"), emit: tre
    tuple val(meta), path("*.rep"), emit: report
    tuple val(meta), path("*.lca"), emit: lca           , optional: true
    tuple val(meta), path("*.all"), emit: all           , optional: true
    tuple val(meta), path("*.unc"), emit: unc           , optional: true
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = meta.single_end ? "--single-reads ${fastqs}" : "--paired-reads ${fastqs}"
    """
    dbprefix=\$(find -L . -name '*.ibf' | sed 's/\\.ibf\$//')

    ganon \\
        classify \\
        --db-prefix \${dbprefix%%.ibf} \\
        $args \\
        --threads $task.cpus \\
        --output-prefix ${prefix} \\
        $input \
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = meta.single_end ? "--single-reads ${fastqs}" : "--paired-reads ${fastqs}"
    """
    touch ${prefix}.tre
    touch ${prefix}.report
    touch ${prefix}.lca
    touch ${prefix}.all
    touch ${prefix}.unc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """
}
