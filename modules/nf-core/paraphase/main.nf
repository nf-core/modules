process PARAPHASE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2f7e7ef6dd6def56bef5252c28a65f279d5fcef09a832510d0de1d6dab09194/data'
        : 'community.wave.seqera.io/library/minimap2_paraphase_samtools:2c52f03fe994efa6'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(config)

    output:
    tuple val(meta), path("*.paraphase.json"), emit: json
    tuple val(meta), path("*.paraphase.bam"), emit: bam
    tuple val(meta), path("*.paraphase.bam.bai"), emit: bai
    tuple val(meta), path("${prefix}_paraphase_vcfs/*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}_paraphase_vcfs/*.vcf.gz.{csi,tbi}"), emit: vcf_index, optional: true
    tuple val("${task.process}"), val('minimap2'), eval("minimap2 --version"), topic: versions, emit: versions_minimap2
    tuple val("${task.process}"), val('paraphase'), eval("paraphase --version"), topic: versions, emit: versions_paraphase
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def config_file = config ? "--config ${config}" : ""
    """
    paraphase \\
        ${args} \\
        --threads ${task.cpus} \\
        --bam ${bam} \\
        --reference ${fasta} \\
        --prefix ${prefix} \\
        ${config_file} \\
        --out .

    for vcf in ${prefix}_paraphase_vcfs/*.vcf;
    do
        bgzip \\
            ${args2} \\
            --threads ${task.cpus} \\
            \$vcf;
        tabix \\
            ${args3} \\
            --threads ${task.cpus} \\
            \$vcf.gz;
    done
    """

    stub:
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def index = args3.contains('--csi') ? 'csi' : 'tbi'
    """
    mkdir ${prefix}_paraphase_vcfs

    touch ${prefix}.paraphase.json
    touch ${prefix}.paraphase.bam
    touch ${prefix}.paraphase.bam.bai
    echo '' | gzip > ${prefix}_paraphase_vcfs/${prefix}_stub.vcf.gz
    touch ${prefix}_paraphase_vcfs/${prefix}_stub.vcf.gz.${index}
    """
}
