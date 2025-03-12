process BBMAP_REPAIR {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'wave.seqera.io/wt/5ce851926d93/wave/build:bbmap-39.18_pigz-2.8--4a8254d490b5baa0' }"

    input:
    tuple val(meta), path(reads)
    val(interleave)

    output:
    tuple val(meta), path("*_repaired.fastq.gz")         , emit: repaired
    tuple val(meta), path("${prefix}_singleton.fastq.gz"), emit: singleton
    path  "versions.yml"                                 , emit: versions
    path  "*.log"                                        , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    in_reads  = ( interleave )  ? "in=${reads[0]}" : "in=${reads[0]} in2=${reads[1]}"
    out_reads = ( interleave )  ? "out=${prefix}_repaired.fastq.gz outs=${prefix}_singleton.fastq.gz"
                                : "out=${prefix}_1_repaired.fastq.gz out2=${prefix}_2_repaired.fastq.gz outs=${prefix}_singleton.fastq.gz"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    repair.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        threads=${task.cpus}
        ${args} \\
        &> ${prefix}.repair.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_1_repaired.fastq.gz
    echo "" | gzip > ${prefix}_2_repaired.fastq.gz
    echo "" | gzip > ${prefix}_singleton.fastq.gz
    touch ${prefix}.repair.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
