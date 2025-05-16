process NEXTDENOVO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa1c1e961de38d24cf36c424a8f4a9920ddd07b63fdb4cfa51c9e3a593c3c979/data' :
        'community.wave.seqera.io/library/flye:2.9.5--d577924c8416ccd8' }"

    input:
    tuple val(meta), path(reads)
    path config

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.stat")     , emit: stat
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${reads} > input.fofn
    nextDenovo \\
        $config \\
        input.fofn \\

    gzip -c ./03.ctg_graph/nd.asm.fasta > ${prefix}.assembly.fasta.gz
  
    mv ./03.ctg_graph/nd.asm.fasta.stat ${prefix}.assembly_info.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$( nextDenovo --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub | gzip -c > ${prefix}.assembly.fasta.gz
    echo contig_1 > ${prefix}.assembly_info.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$( nextDenovo --version )
    END_VERSIONS
    """
}
