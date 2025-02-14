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
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = aln.name.endsWith(".gz")
    def aln_name = aln.name
    if (is_compressed) {
        aln_name = aln.name.replace(".gz", "")
    }
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $aln > $aln_name
    fi

    reformat.pl \\
        $args \\
        ${informat} \\
        ${outformat} \\
        ${aln_name} \\
        ${prefix}.${outformat}

    gzip ${prefix}.${outformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hhsuite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.${outformat}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hhsuite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
