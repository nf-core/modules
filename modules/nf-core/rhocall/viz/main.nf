process RHOCALL_VIZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7dbf7021085cfea72a20cafffe57fcf47392706d9a433f1143f1e60b389b85ae/data':
        'community.wave.seqera.io/library/rhocall:0.5.1--a7eced77e39d2b82' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(roh)

    output:
    tuple val(meta), path("${prefix}/${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}/${prefix}.wig"), emit: wig
    tuple val("${task.process}"), val("rhocall"), eval("rhocall --version | sed 's/rhocall, version //'"), topic: versions, emit: versions_rhocall

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rhocall \\
        viz \\
        $args \\
        -r $roh \\
        --out_dir ${prefix} \\
        $vcf

    mv ${prefix}/output.bed ${prefix}/${prefix}.bed
    mv ${prefix}/output.wig ${prefix}/${prefix}.wig
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.bed
    touch ${prefix}/${prefix}.wig
    """
}
