process HHSUITE_REFORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hhsuite:3.3.0--py39pl5321h0dd7abe_13':
        'biocontainers/hhsuite:3.3.0--py39pl5321h0dd7abe_13' }"

    input:
    tuple val(meta), path(aln)
    val(informat)
    val(outformat)

    output:
    tuple val(meta), path("${prefix}.${outformat}.gz"), emit: msa
    tuple val("${task.process}"), val("hhsuite"), eval("hhblits -h 2>&1 | sed -n '1s/^HHblits \\([0-9.]\\+\\):/\\1/p'"), topic: versions, emit: versions_hhsuite

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = aln.getExtension() == "gz" ? true : false
    def aln_name = is_compressed ? aln.getBaseName() : aln
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${aln} > ${aln_name}
    fi

    reformat.pl \\
        $args \\
        ${informat} \\
        ${outformat} \\
        ${aln_name} \\
        ${prefix}.${outformat}

    gzip ${prefix}.${outformat}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo "" | gzip > ${prefix}.${outformat}.gz
    """
}
