process SPACERANGER_MKREF {
    tag "$fasta"
    label 'process_high'

    // TODO push to nf-core docker
    container "ghcr.io/grst/spaceranger:2.1.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "SPACERANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path fasta
    path gtf
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    spaceranger \\
        mkref \\
        --genome=$reference_name \\
        --fasta=$fasta \\
        --genes=$gtf \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: \$(spaceranger -V | sed -e "s/spaceranger spaceranger-//g")
    END_VERSIONS
    """
}
