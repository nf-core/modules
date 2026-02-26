process HMMCOPY_READCOUNTER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b656ae8c47699d1ae8f4d97995da8fda08aaed7a23d9b7078e034ff0981f7487/data' :
        'community.wave.seqera.io/library/hmmcopy_samtools:875db3767c6d4ea2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.wig"), emit: wig
    tuple val("${task.process}"), val('hmmcopy'), eval("echo 0.1.1"), topic: versions, emit: versions_hmmcopy
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    // Note that piping the cram directly into the tool didn't work.
    def convert_cram     = bam.Extension == "cram" ? "samtools view -T ${fasta} -h ${bam} -o temp.bam##idx##temp.bam.bai --write-index " : ""
    def input            = bam.Extension == "cram" ? "temp.bam" : "${bam}"
    def cleanup          = bam.Extension == "cram" ? "rm temp.bam{,.bai}" : ""
    """
    ${convert_cram}

    readCounter \\
        $args \\
        ${input} > ${prefix}.wig

    ${cleanup}
    """
    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.wig
    """
}
