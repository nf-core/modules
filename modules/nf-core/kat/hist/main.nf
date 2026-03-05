process KAT_HIST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kat:2.4.2--py38hfc5f9d8_2':
        'biocontainers/kat:2.4.2--py38hfc5f9d8_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.hist")                   , emit: hist
    tuple val(meta), path("*.hist.dist_analysis.json"), emit: json
    tuple val(meta), path("*.png")                    , emit: png           , optional: true
    tuple val(meta), path("*.ps")                     , emit: ps            , optional: true
    tuple val(meta), path("*.pdf")                    , emit: pdf           , optional: true
    tuple val(meta), path("*-hash.jf*")               , emit: jellyfish_hash, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/merqury/merqury

Reason:
This module no longer works in conda due to glibc incompatibilities with plotting libraries
This module is no longer maintained by the authors
"""
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert false: deprecation_message
    """
    kat hist \\
        --threads $task.cpus \\
        --output_prefix ${prefix}.hist \\
        $args \\
        $reads

    ls -l

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kat: \$( kat hist --version | sed 's/kat //' )
    END_VERSIONS
    """

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/merqury/merqury

Reason:
This module no longer works in conda due to glibc incompatibilities with plotting libraries
This module is no longer maintained by the authors
"""
    def prefix    = task.ext.prefix ?: "${meta.id}"
    assert false: deprecation_message
    """
    touch ${prefix}.hist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kat: \$( kat hist --version | sed 's/kat //' )
    END_VERSIONS
    """
}
