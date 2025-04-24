process HMMER_HMMBUILD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07c4cbd91c4459dc86b13b5cd799cacba96b27d66c276485550d299c7a4c6f8a/data' :
        'community.wave.seqera.io/library/hmmer:3.4--cb5d2dd2e85974ca' }"

    input:
    tuple val(meta), path(alignment)
    path mxfile

    output:
    tuple val(meta), path("*.hmm.gz"), emit: hmm
    path "*.hmmbuild.txt",             emit: hmmbuildout
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def mxfileopt = mxfile ? "--mxfile ${mxfile}" : ""

    """
    hmmbuild \\
        $args \\
        --cpu $task.cpus \\
        -n ${prefix}  \\
        -o ${prefix}.hmmbuild.txt \\
        ${mxfileopt} \\
        ${prefix}.hmm \\
        $alignment

    gzip ${prefix}.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(hmmbuild -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.hmm.gz
    touch ${prefix}.hmmbuild.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(hmmbuild -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
    END_VERSIONS
    """
}
