process MUDSKIPPER_BULK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mudskipper:0.1.0--h9f5acd7_1':
        'biocontainers/mudskipper:0.1.0--h9f5acd7_1' }"

    input:
    tuple val(meta), path(bam)
    path index
    path gtf
    val rad

    output:
    tuple val(meta), path("${prefix}.bam"), optional:true, emit: bam
    tuple val(meta), path("${prefix}.rad"), optional:true, emit: rad
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.transcriptome"
    def annot_param = ""
    if (index) {
        annot_param = "--index ${index}"
    } else {
        annot_param = "--gtf ${gtf}"
    }
    def suffix = rad ? "rad" : "bam"
    """
    export RUST_BACKTRACE=full
    mudskipper \\
        bulk \\
        ${annot_param} \\
        --alignment ${bam} \\
        --out ${prefix}.${suffix} \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mudskipper: \$(echo \$(mudskipper -V 2>&1) | sed 's/^.*mudskipper //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.transcriptome"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mudskipper: \$(echo \$(mudskipper -V 2>&1) | sed 's/^.*mudskipper //' )
    END_VERSIONS
    """
}
