process ANGSD_NGSADMIX {
    tag "${meta.id}_K${k}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--h13024bc_4':
        'quay.io/biocontainers/angsd:0.940--h13024bc_4' }"

    input:
    tuple val(meta), path(beagle_file), val(k)

    output:
    tuple val(meta), path("*.qopt")   , emit: qopt
    tuple val(meta), path("*.fopt.gz"), emit: fopt
    tuple val(meta), path("*.log")    , emit: log
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.K${k}"

    """
    NGSadmix \\
        -likes ${beagle_file} \\
        -K ${k} \\
        -P ${task.cpus} \\
        -outfiles ${prefix} \\
        ${args}
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.K${k}"
    """
    echo ${args}

    touch ${prefix}.qopt
    echo "" | gzip > ${prefix}.fopt.gz
    touch ${prefix}.log
    """
}
