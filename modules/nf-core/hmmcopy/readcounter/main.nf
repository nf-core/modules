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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    // Note that piping the cram directly into the tool didn't work.
    def convert_cram     = bam.Extension == "cram" ? "samtools view -T ${fasta} -h ${bam} -o temp.bam##idx##temp.bam.bai --write-index " : ""
    def input            = bam.Extension == "cram" ? "temp.bam" : "${bam}"
    def cleanup          = bam.Extension == "cram" ? "rm temp.bam{,.bai}" : ""
    def samtools_version = bam.Extension == "cram" ? "samtools: \$(samtools --version |& sed '1!d ; s/samtools //')" : ""
    """
    ${convert_cram}

    readCounter \\
        $args \\
        ${input} > ${prefix}.wig

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: $VERSION
        $samtools_version
    END_VERSIONS
    """
    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def samtools_version = bam.Extension == "cram" ? "samtools: \$(samtools --version |& sed '1!d ; s/samtools //')" : ""
    """
    touch ${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: $VERSION
        $samtools_version
    END_VERSIONS
    """
}
