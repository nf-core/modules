process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'
    ext prefix: "${meta.id}", args: ''

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam, stageAs: "?/*")
    tuple val(meta2), path(fasta), path(fai)
    val index_format

    output:
    tuple val(meta), path("${task.ext.prefix}.bam"), emit: bam, optional: true
    tuple val(meta), path("${task.ext.prefix}.cram"), emit: cram, optional: true
    tuple val(meta), path("${task.ext.prefix}.sam"), emit: sam, optional: true
    tuple val(meta), path("${task.ext.prefix}.${extension}.{crai,csi,bai}"), emit: index, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when: 
    task.ext.when == null || task.ext.when

    script:
    extension = task.ext.args.contains("--output-fmt sam")
        ? "sam"
        : task.ext.args.contains("--output-fmt cram")
            ? "cram"
            : "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    //setting default values
    def write_index = ""
    def output_file = "${task.ext.prefix}.${extension}"

    // Update if index is requested
    if (index_format != '' && index_format) {
        write_index = "--write-index"
        output_file = "${task.ext.prefix}.${extension}##idx##${task.ext.prefix}.${extension}.${index_format}"
    }
    def is_sam = (bam instanceof List ? bam[0] : bam).name.endsWith('.sam')
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (extension == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }
    if ("${bam}" == "${task.ext.prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    if ("${bam}" == "${task.ext.prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    def input_source = is_sam ? "${bam}" : "-"
    def pre_command = is_sam ? "" : "samtools cat ${bam} | "

    """
    ${pre_command}samtools sort \\
        ${task.ext.args} \\
        -T ${task.ext.prefix} \\
        --threads ${task.cpus} \\
        ${reference} \\
        -o ${output_file} \\
        ${write_index} \\
        ${input_source}
    """

    stub:
    extension = task.ext.args.contains("--output-fmt sam")
        ? "sam"
        : task.ext.args.contains("--output-fmt cram")
            ? "cram"
            : "bam"

    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (extension == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }

    index = index_format ? "touch ${task.ext.prefix}.${extension}.${index_format}" : ""

    """
    touch ${task.ext.prefix}.${extension}
    ${index}
    """
}
