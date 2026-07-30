process ANGSD_SFS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--h13024bc_4':
        'quay.io/biocontainers/angsd:0.940--h13024bc_4' }"

    input:    
    tuple val(meta), path(pop1_saf_idx), path(pop1_saf_pos), path(pop1_saf)
    tuple val(meta2), path(pop2_saf_idx), path(pop2_saf_pos), path(pop2_saf) // Optional: use for 2-dimensional SFS estimation

    output:
    tuple val(meta_out), path("*.sfs"), emit: sfs
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_2d = meta2 as boolean
    def meta_out = is_2d ? [id: "${meta.id}_${meta2.id}", pop1: meta.pop, pop2: meta2.pop] : meta
    def pop2_input = is_2d ? "${pop2_saf_idx}" : ''
    def prefix = task.ext.prefix ?: meta_out.id
    
    """
    realSFS \\
        ${pop1_saf_idx} \\
        ${pop2_input} \\
        ${args} \\
        > ${prefix}.sfs
    """

    stub:
    def prefix = task.ext.prefix ?: meta_out.id
    """
    touch ${prefix}.sfs
    """
}
