process SCRAMBLE_CLUSTERIDENTIFIER {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::scramble=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scramble:1.0.1--h779adbc_1':
        'quay.io/biocontainers/scramble:1.0.1--h779adbc_1' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path fasta

    output:
    tuple val(meta), path("*.clusters.txt") , emit: clusters
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
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
}
