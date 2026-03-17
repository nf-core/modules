process CELLRANGER_MKREF {
    tag "$fasta"
    label 'process_high'

    container "nf-core/cellranger:10.0.0"

    input:
    path fasta
    path gtf
    val reference_name

    output:
    path "${reference_name}", emit: reference
    tuple val("${task.process}"), val('cellranger'), eval('cellranger --version | sed "s/.*-//"'), emit: versions_cellranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    // --localcores is passed to the martian runtime and specifies the number of allocated jobs
    // --nthreads is passed to the STAR index generation.
    // see also https://github.com/nf-core/scrnaseq/issues/329
    """
    cellranger \\
        mkref \\
        --genome=$reference_name \\
        --fasta=$fasta \\
        --genes=$gtf \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --nthreads=${task.cpus} \\
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir $reference_name
    touch ${reference_name}/empty_file
    """

}
