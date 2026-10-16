process VARDICTRUST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/12/12806f94133378ff8728b41f16ce19a669f61750f18db4ba3563af807ad45cc6/data'
:         'community.wave.seqera.io/library/htslib_vardictrust:913eef801f19dc83' }"

    input:
    tuple val(meta), path(bams), path(bais), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('vardictrust'), eval('vardictrust --version | head -n 1 | sed \'s/.* //\''), emit: versions_vardictrust, topic: versions
    tuple val("${task.process}"), val('htslib'), eval("tabix -h 2>&1 | sed -n '2s/Version: *//p'"), emit: versions_htslib, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-c 1 -S 2 -E 3 -z'
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // A lone BAM may arrive as a 1-element List; unwrap it so it never mis-dispatches as somatic
    def bam_input = bams instanceof List && bams.size() == 1 ? bams[0] : bams
    def somatic = bam_input instanceof List && bam_input.size() == 2
    // Single-quoted Groovy: the emitted shell keeps -b "tumor|normal" double-quoted so | is never a pipe
    def input = somatic ? '-b "' + bam_input[0] + '|' + bam_input[1] + '"' : "-b ${bam_input}"
    """
    set -o pipefail

    vardictrust \\
        ${args} \\
        ${args2} \\
        ${input} \\
        -N "${meta.id}" \\
        --threads ${task.cpus} \\
        -G ${fasta} \\
        ${bed} \\
        | bgzip ${args3} --threads ${task.cpus} > ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo '' | gzip > ${prefix}.vcf.gz
    """
}
