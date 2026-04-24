process MODKIT_ENTROPY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    // stageAs '?/*' prevents filename collisions when multiple BAMs from the
    // same sample (e.g. technical replicates) are passed to a single run.
    tuple val(meta),  path(bams, stageAs: "in/?/*"), path(bais, stageAs: "in/?/*")
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(regions)

    output:
    tuple val(meta), path("*.bed")                                  , emit: bed           , optional: true
    tuple val(meta), path("entropy_regions/*.bed")                  , emit: regions_bed   , optional: true
    tuple val(meta), path("entropy_regions/*.bedgraph")             , emit: bedgraph      , optional: true
    tuple val(meta), path("entropy_regions/*.tsv")                  , emit: tsv           , optional: true
    tuple val(meta), path("*.log")                                  , emit: log           , optional: true
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def bam_args = bams instanceof List ? bams.collect { "--in-bam ${it}" }.join(' ') : "--in-bam ${bams}"
    // modkit entropy's --out-bed expects a FILE without --regions, and a DIRECTORY with --regions
    def out_arg  = regions ? "--regions ${regions} --out-bed entropy_regions --prefix ${prefix}" : "--out-bed ${prefix}.bed"
    def mkdir    = regions ? "mkdir -p entropy_regions" : ""
    """
    ${mkdir}

    modkit \\
        entropy \\
        $args \\
        --threads ${task.cpus} \\
        --ref ${fasta} \\
        ${out_arg} \\
        ${bam_args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
