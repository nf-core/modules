process RHOCALL_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7dbf7021085cfea72a20cafffe57fcf47392706d9a433f1143f1e60b389b85ae/data':
        'community.wave.seqera.io/library/rhocall:0.5.1--a7eced77e39d2b82' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(roh)
    path bed

    output:
    tuple val(meta), path("*_rhocall.vcf"), emit: vcf
    tuple val("${task.process}"), val("rhocall"), eval("rhocall --version | sed 's/rhocall, version //'"), topic: versions, emit: versions_rhocall

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def az_bed = bed ? "-b ${bed}" : ''
    """
    export MPLCONFIGDIR=\$PWD

    rhocall \\
        annotate \\
        $args \\
        $az_bed \\
        -r $roh \\
        -o ${prefix}_rhocall.vcf \\
        $vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD

    touch ${prefix}_rhocall.vcf
    """
}
