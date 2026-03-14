process ANGSD_CONTAMINATION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--hf5e1c6e_3':
        'biocontainers/angsd:0.940--hf5e1c6e_3' }"

    input:
    tuple val(meta), path(icounts)
    tuple val(meta2), path(hapmap_file)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def seed_cmd = args.contains("-s ") ? '' : '-s 1'
    """
    contamination \
        ${args} \
        ${seed_cmd} \
        -a ${icounts} \
        -h ${hapmap_file} \
        -p ${task.cpus} \
        2>| >(tee ${prefix}.txt >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
