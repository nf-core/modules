process CELLRANGER_MKREF {
    tag "$fasta"
    label 'process_high'

    container "nf-core/cellranger:8.0.0"

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    cellranger \\
        mkref \\
        --genome=$reference_name \\
        --fasta=$fasta \\
        --genes=$gtf \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    mkdir $reference_name
    touch ${reference_name}/empty_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

}
