process CELLRANGER_MKVDJREF {
    tag "$fasta"
    label 'process_high'

    container "nf-core/cellranger:8.0.0"

    input:
    path fasta          // optional
    path gtf            // optional
    path seqs           // optional
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKVDJREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args        = task.ext.args ?: ''
    def gtf_in      = gtf           ? "--genes ${gtf}"      : ""
    def fasta_in    = fasta         ? "--fasta ${fasta}"    : ""
    def seqs_in     = seqs          ? "--seqs ${seqs}"      : ""

    """
    cellranger \\
        mkvdjref \\
        --genome=$reference_name \\
        ${gtf_in} \\
        ${fasta_in} \\
        ${seqs_in} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir ${reference_name}
    mkdir ${reference_name}/fasta
    echo stub > ${reference_name}/fasta/regions.fa
    echo stub > ${reference_name}/reference.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
