process RASTAIR_CALL {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15120636da858ba73a2493281bfa418005f08c0ed09369a837c05f3f9e14a4a6/data' :
        'community.wave.seqera.io/library/rastair:0.8.2--bf70eeab4121509c' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), val(parsed_trim_OT)
    tuple val(meta6), val(parsed_trim_OB)

    output:
    tuple val(meta), path("*.rastair_call.txt"),    emit: txt
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nt_OT_to_trim = meta.trim_OT ?: parsed_trim_OT
    def nt_OB_to_trim = meta.trim_OB ?: parsed_trim_OB

    """
    rastair call \\
        --threads ${task.cpus} \\
        --nOT ${nt_OT_to_trim} \\
        --nOB ${nt_OB_to_trim} \\
        --fasta-file ${fasta} \\
        ${bam} > ${prefix}.rastair_call.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_call.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version 2>&1 || echo "stub")
    END_VERSIONS
    """
}
