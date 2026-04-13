process TELOMEREHUNTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80328294e56cd32cb354e132be9fe29b20e59acbdd1e071cd94aec5c21f9abda/data'
        : 'community.wave.seqera.io/library/python_pysam_samtools_numpy_pruned:e5e0b9c3eb477e2e' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(control_bam), path(control_bai)
    tuple val(meta2), path(fasta), path(fai), path(cytoband)

    output:
    tuple val(meta), path("${prefix}/${prefix}_summary.tsv")        , emit: summary
    tuple val(meta), path("${prefix}/tumor_TelomerCnt_${prefix}/")  , emit: tumor
    tuple val(meta), path("${prefix}/control_TelomerCnt_${prefix}/"), emit: control, optional: true
    tuple val("${task.process}"), val('telomerehunter'), eval("pip show telomerehunter | sed -n 's/^Version: //p'"), emit: versions_telomerehunter, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version |& sed '1!d; s/^.*samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cytoband_arg = cytoband ? "-b ${cytoband}" : ""
    // telomerehunter doesn't support CRAM (pysam opened in BAM-only mode)
    def tumor_is_cram = tumor_bam.name.endsWith(".cram")
    def control_is_cram = control_bam ? control_bam.name.endsWith(".cram") : false
    def tumor = tumor_is_cram ? "tumor.bam" : "${tumor_bam}"
    def control = control_bam ? (control_is_cram ? "control.bam" : "${control_bam}") : ""
    def tumor_to_bam_cmd = tumor_is_cram ? "samtools view -T ${fasta} -b -o tumor.bam ${tumor_bam} && samtools index tumor.bam" : ""
    def control_to_bam_cmd = control_is_cram ? "samtools view -T ${fasta} -b -o control.bam ${control_bam} && samtools index control.bam" : ""
    """
    ${tumor_to_bam_cmd}
    ${control_to_bam_cmd}

    telomerehunter \\
        -ibt ${tumor} \\
        ${control ? "-ibc ${control}" : ""} \\
        ${cytoband_arg} \\
        -o . \\
        -p ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def control_stub = control_bam ? """
    mkdir -p ${prefix}/control_TelomerCnt_${prefix}
    touch ${prefix}/control_TelomerCnt_${prefix}/${prefix}_control_summary.tsv
    """ : ""
    """
    mkdir -p ${prefix}/tumor_TelomerCnt_${prefix}
    touch ${prefix}/${prefix}_summary.tsv
    touch ${prefix}/tumor_TelomerCnt_${prefix}/${prefix}_tumor_summary.tsv
    ${control_stub}
    """
}
