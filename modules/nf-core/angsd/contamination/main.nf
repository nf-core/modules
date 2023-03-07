process ANGSD_CONTAMINATION {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::angsd=0.940"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'quay.io/biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(icounts)
    tuple val(meta2), path(hapmap_file)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    contamination \
        ${args} \
        -a ${icounts} \
        -h ${hapmap_file} \
        -p ${task.cpus} \
        2> ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}
