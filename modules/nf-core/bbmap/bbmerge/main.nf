process BBMAP_BBMERGE {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    val(interleave)

    output:
    tuple val(meta), path("*_merged.fastq.gz")  , emit: merged
    tuple val(meta), path("*_unmerged.fastq.gz"), emit: unmerged
    tuple val(meta), path("*_ihist.txt")        , emit: ihist
    path  "versions.yml"                        , emit: versions
    path  "*.log"                               , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_reads = ( interleave ) ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def out_reads = ( interleave ) ? "out=${prefix}_merged.fastq.gz outu=${prefix}_unmerged.fastq.gz" : "out=${prefix}_merged.fastq.gz outu1=${prefix}_1_unmerged.fastq.gz outu2=${prefix}_2_unmerged.fastq.gz"

    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbmerge.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        ihist=${prefix}_ihist.txt \\
        $args \\
        &> ${prefix}.bbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_files = ( interleave ) ? "${prefix}_merged.fastq.gz ${prefix}_unmerged.fastq.gz" : "${prefix}_merged.fastq.gz ${prefix}_1_unmerged.fastq.gz ${prefix}_2_unmerged.fastq.gz"

    """
    echo "" | gzip | tee $out_files
    touch ${prefix}_ihist.txt
    touch ${prefix}.bbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
