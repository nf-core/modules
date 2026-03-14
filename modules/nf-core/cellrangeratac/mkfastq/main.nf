process CELLRANGERATAC_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    container "nf-core/cellranger-atac-mkfastq:2.1.0"

    input:
    path bcl
    path csv

    output:
    path "${bcl.getSimpleName()}/outs/fastq_path/*.fastq.gz"  , emit: fastq
    tuple val("${task.process}"), val('cellrangeratac'), eval("cellranger-atac --version | sed 's/.*cellranger-atac-//'"), emit: versions_cellrangeratac, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGERATAC_MKFASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    """
    cellranger-atac mkfastq --id=${bcl.getSimpleName()} \
        --run=$bcl \
        --csv=$csv \
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGERATAC_MKFASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    """
    mkdir -p "${bcl.getSimpleName()}/outs/fastq_path/"
    echo "" | gzip > ${bcl.getSimpleName()}/outs/fastq_path/fake_file.fastq.gz
    """
}
