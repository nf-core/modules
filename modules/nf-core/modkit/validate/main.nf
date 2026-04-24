process MODKIT_VALIDATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    // Multiple sample pairs can be supplied as list inputs; BAMs and BEDs are
    // paired by index and passed as repeated `--bam-and-bed <BAM> <BED>`.
    tuple val(meta),  path(bams, stageAs: "bams/?/*"), path(bais, stageAs: "bais/?/*"), path(truth_beds, stageAs: "truth/?/*")

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def bam_list   = bams       instanceof List ? bams       : [bams]
    def bed_list   = truth_beds instanceof List ? truth_beds : [truth_beds]
    if (bam_list.size() != bed_list.size()) {
        error "MODKIT_VALIDATE: bams and truth_beds must have the same length (got ${bam_list.size()} vs ${bed_list.size()})"
    }
    def pair_args  = [bam_list, bed_list].transpose().collect { b, g -> "--bam-and-bed ${b} ${g}" }.join(' ')
    """
    modkit \\
        validate \\
        $args \\
        --threads ${task.cpus} \\
        ${pair_args} \\
        --out-filepath ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
