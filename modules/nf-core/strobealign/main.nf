process STROBEALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59cce6872df48a1e5cc9ccee89f066210694c6ec9f62d9c931cc6925ca0f6a5f/data' :
        'community.wave.seqera.io/library/htslib_samtools_strobealign_pigz:4fa4f439c6bea386' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    val sort_bam

    output:
    tuple val(meta), path("*.bam")      , emit: bam , optional: true
    tuple val(meta), path("*.cram")     , emit: cram, optional: true
    tuple val(meta), path("*.csi")      , emit: csi , optional: true
    tuple val(meta), path("*.crai")     , emit: crai, optional: true
    tuple val(meta), path("*.paf.gz")   , emit: paf , optional: true
    tuple val(meta), path("*.tsv.gz")   , emit: tsv , optional: true
    tuple val(meta), path("*.sti")      , emit: sti , optional: true
    tuple val("${task.process}"), val('strobealign'), eval("strobealign --version"), topic: versions, emit: versions_strobealign
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //'"), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Determine output file extension and command
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args.contains("-x")
        ? "paf"
        : args.contains("--aemb")
            ? "tsv"
            : args.contains("--create-index")
                ? "sti"
                : args2.contains("--output-fmt cram")
                    ? "cram"
                    : sort_bam && args2.contains("-O cram")
                        ? "cram"
                        : !sort_bam && args2.contains("-C")
                            ? "cram"
                            : "bam"
    def reference = fasta && extension == "cram" ? "--reference ${fasta}" : ""
    if (!fasta && extension == "cram") {
        error("Fasta reference is required for CRAM output")
    }
    def output_cmd = extension == "bam" || extension == "cram"
        ? "samtools ${samtools_command} ${args2} ${reference} --threads ${task.cpus} -o ${prefix}.${extension} -"
        : extension == "paf"
            ? "pigz ${args2} > ${prefix}.paf.gz"
            : extension == "tsv"
                ? "pigz ${args2} > ${prefix}.tsv.gz"
                : extension == "sti"
                    ? "tee ${prefix}.log > /dev/null"
                    : error("Unable to determine output command for extension: ${extension}")

    """
    strobealign \\
        ${args} \\
        -t ${task.cpus} \\
        ${fasta} \\
        ${reads} \\
        | ${output_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-x")
        ? "paf"
        : args.contains("--aemb")
            ? "tsv"
            : args2.contains("--output-fmt cram")
                ? "cram"
                : sort_bam && args2.contains("-O cram")
                    ? "cram"
                    : !sort_bam && args2.contains("-C")
                        ? "cram"
                        : "bam"

    """
    touch ${prefix}.${extension}
    touch ${prefix}.csi
    touch ${prefix}.crai
    echo "" | pigz > ${prefix}.paf.gz
    echo "" | pigz > ${prefix}.tsv.gz
    touch ${prefix}.sti
    """
}
