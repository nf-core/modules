process VG_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.xg")       , emit: xg
    tuple val(meta), path("*.vgi")      , emit: vg_index, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vg index \\
        --temp-dir . \\
        --threads ${task.cpus} \\
        --xg-name ${prefix}.xg \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def vg_index = args.contains('--index-sorted-vg') ? "touch ${prefix}.vg.vgi" : ""

    """
    touch ${prefix}.xg
    ${vg_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """
}
