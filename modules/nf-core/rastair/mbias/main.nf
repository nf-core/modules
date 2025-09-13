process RASTAIR_MBIAS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f0a47f3c0c4f521ed0623cd709c51fb1ece4df1fb4bd85c75d04e0383a8c5d4/data' :
        'community.wave.seqera.io/library/rastair:0.8.2--b09a8e25a0d53059' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.rastair_mbias.txt"),   emit: txt
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rastair mbias \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        ${bam} > ${prefix}.rastair_mbias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
