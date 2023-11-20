process PARABRICKS_INDEXGVCF {
    tag "$meta.id"
    label 'process_high'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"

    /*
    NOTE: Parabricks requires the files to be non-symlinked
    Do not change the stageInMode to soft linked! This is default on Nextflow.
    If you change this setting be careful.
    */
    stageInMode "copy"

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("*.vcf.tbi"), emit: gvcf_index
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pbrun \\
        indexgvcf \\
        --input $gvcf \\
        --num-threads $task.cpus \\
        --num-gpus $task.accellerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.g.vcf.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
