process MODKIT_DMR_PAIR {
    tag "${meta.id}_vs_${meta2.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    tuple val(meta),  path(control_bed,    stageAs: 'control/*'), path(control_tbi,    stageAs: 'control/*')
    tuple val(meta2), path(experiment_bed, stageAs: 'experiment/*'), path(experiment_tbi, stageAs: 'experiment/*')
    tuple val(meta3), path(fasta),          path(fai)
    tuple val(meta4), path(regions_bed)

    output:
    tuple val(meta), path("*.dmr_sites.bed")   , emit: sites
    tuple val(meta), path("*.dmr_segments.bed"), emit: segments, optional: true
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def regions_arg = regions_bed     ? "--regions-bed ${regions_bed}" : ''
    def segment_arg = "--segment ${prefix}.dmr_segments.bed"
    """
    modkit \\
        dmr pair \\
        -a ${control_bed} \\
        -b ${experiment_bed} \\
        --ref ${fasta} \\
        --base C \\
        --header \\
        --threads ${task.cpus} \\
        ${regions_arg} \\
        ${segment_arg} \\
        ${args} \\
        --out-path ${prefix}.dmr_sites.bed \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dmr_sites.bed
    touch ${prefix}.dmr_segments.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """
}
