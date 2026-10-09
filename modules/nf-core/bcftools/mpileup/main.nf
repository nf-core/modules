process BCFTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(bam), path(intervals_mpileup, stageAs: 'mpileup_intervals/*'), path(intervals_call, stageAs: 'call_intervals/*')
    tuple val(meta2), path(fasta), path(fai)
    val save_mpileup

    output:
    tuple val(meta), path("*.{bcf,vcf}{,.gz}"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}"), emit: index, optional: true
    tuple val(meta), path("*stats.txt"), emit: stats
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def bgzip_mpileup = save_mpileup ? "bgzip ${prefix}.mpileup" : ""
    def intervals_mpileup_cmd = intervals_mpileup ? "-T ${intervals_mpileup}" : ""
    def intervals_call_cmd = intervals_call ? "-T ${intervals_call}" : ""
    def extension = args3.contains("--output-type b") || args3.contains("-Ob")
        ? "bcf.gz"
        : args3.contains("--output-type u") || args3.contains("-Ou")
            ? "bcf"
            : args3.contains("--output-type z") || args3.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args3.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    """
    echo "${meta.id}" > sample_name.list

    bcftools \\
        mpileup \\
        --fasta-ref ${fasta} \\
        ${args} \\
        ${bam} \\
        ${intervals_mpileup_cmd} \\
        ${mpileup} \\
        | bcftools call --output-type v ${args2} ${intervals_call_cmd} \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${prefix}.${extension} ${args3}

    ${bgzip_mpileup}

    bcftools stats ${prefix}.${extension} > ${prefix}.bcftools_stats.txt
    """

    stub:
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args3.contains("--output-type b") || args3.contains("-Ob")
        ? "bcf.gz"
        : args3.contains("--output-type u") || args3.contains("-Ou")
            ? "bcf"
            : args3.contains("--output-type z") || args3.contains("-Oz")
                ? "vcf.gz"
                : args3.contains("--output-type v") || args3.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index = args3.contains("--write-index=tbi") || args3.contains("-W=tbi")
        ? "tbi"
        : args3.contains("--write-index=csi") || args3.contains("-W=csi")
            ? "csi"
            : args3.contains("--write-index") || args3.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""
    def create_mpileup = save_mpileup ? "echo '' | gzip > ${prefix}.mpileup.gz" : ""
    """
    touch ${prefix}.bcftools_stats.txt
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    ${create_mpileup}
    """
}
