process BWAMEM3_MEM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80064f4b7f3a9d932fd4019f39a56ed949d403765b6a11a849c7ce19e74490ed/data'
        : 'community.wave.seqera.io/library/bwa-mem3_samtools_htslib:c28a809633c294ed'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val sort_bam

    output:
    tuple val(meta), path("*.sam"), emit: sam, optional: true
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bwamem3'), eval('bwa-mem3 version 2>/dev/null | head -n1 | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem3, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'

    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher = (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension == "cram" ? "--reference ${fasta}" : ""
    if (!fasta && extension == "cram") {
        error("Fasta reference is required for CRAM output")
    }
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem3 \\
        mem \\
        ${args} \\
        -t ${task.cpus} \\
        \$INDEX \\
        ${reads} \\
        | samtools ${samtools_command} ${args2} -@ ${task.cpus} ${reference} -o ${prefix}.${extension} -
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher = (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    if (!fasta && extension == "cram") {
        error("Fasta reference is required for CRAM output")
    }

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    }
    else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }
    """
    touch ${prefix}.${extension}
    ${create_index}
    """
}
