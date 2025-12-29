process SCRAMBLE_CLUSTERIDENTIFIER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65d3a32dfd347b370e87589189717c75468e6d737b7cee6931e4dae21ce1a9cf/data':
        'community.wave.seqera.io/library/bioconductor-pwalign_scramble:31d27d3832b0689e' }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.clusters.txt") , emit: clusters
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    // The tool does not contain a way to specify the reference file when using CRAM files.
    // It just looks in the header of the CRAM file where the reference file is located,
    // but that reference can't always be fetched since most test data is created on
    // another machine. I had to find another way to specify the reference and I
    // found that I could create an md5 cache of a specified fasta and supply it to
    // the REF_PATH environment variable. This way the tool uses the correct reference.
    // An issue has been made about this: https://github.com/GeneDx/scramble/issues/27
    // The reference code is a placeholder until this issue has been fixed.
    def reference = fasta ? "wget https://raw.githubusercontent.com/samtools/samtools/master/misc/seq_cache_populate.pl && perl seq_cache_populate.pl -root ./md5_ref ${fasta} && export REF_PATH=`pwd`/md5_ref/%2s/%2s/%s" : ""
    """
    ${reference}

    cluster_identifier \\
        ${args} \\
        ${input} \\
        > ${prefix}.clusters.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scramble: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.clusters.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scramble: ${VERSION}
    END_VERSIONS
    """
}
