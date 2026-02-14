process BWAMEM2_INDEX {
    tag "${fasta}"
    // NOTE Requires 28N GB memory where N is the size of the reference sequence, floor of 280M
    // source: https://github.com/bwa-mem2/bwa-mem2/issues/9
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data'
        : 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}"), emit: index
    tuple val("${task.process}"), val('bwamem2'), eval('bwa-mem2 version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id ? "${meta.id}" : "${fasta.baseName}"
    def args = task.ext.args ?: ''
    def suffix = task.ext.suffix ?: "${fasta.baseName}"
    """
    mkdir ${prefix}
    bwa-mem2 \\
        index \\
        ${args} \\
        -p ${prefix}/${suffix} \\
        ${fasta}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id ? "${meta.id}" : "${fasta.baseName}"
    def suffix = task.ext.suffix ?: "${fasta.baseName}"
    """
    mkdir ${prefix}
    touch ${prefix}/${suffix}.0123
    touch ${prefix}/${suffix}.ann
    touch ${prefix}/${suffix}.pac
    touch ${prefix}/${suffix}.amb
    touch ${prefix}/${suffix}.bwt.2bit.64
    """
}
