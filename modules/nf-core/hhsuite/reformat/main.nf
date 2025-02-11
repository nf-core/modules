process HHSUITE_REFORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/60413102dd36ffb8767525d16df897c5b5cee15a8af56bbcf10f407537aa823b/data':
        'community.wave.seqera.io/library/hhsuite_perl:5e6367cac8ba3a53' }"

    input:
    tuple val(meta), path(aln)
    val(informat)
    val(outformat)

    output:
    tuple val(meta), path("${prefix}.${outformat}"), emit: msa
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = aln.name.endsWith(".gz")
    def aln_name = aln.name
    if (is_compressed)
        aln_name = aln.name.replace(".gz", "")
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hhsuite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${outformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hhsuite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
