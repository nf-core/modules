process BWA_ALN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sai"), emit: sai
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), topic: versions, emit: versions_bwa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

        bwa aln \\
            $args \\
            -t $task.cpus \\
            -f ${prefix}.sai \\
            \$INDEX \\
            ${reads}
        """
    } else {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

        bwa aln \\
            $args \\
            -t $task.cpus \\
            -f ${prefix}.1.sai \\
            \$INDEX \\
            ${reads[0]}

        bwa aln \\
            $args \\
            -t $task.cpus \\
            -f ${prefix}.2.sai \\
            \$INDEX \\
            ${reads[1]}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        touch ${prefix}.sai
        """
    } else {
        """
        touch ${prefix}.1.sai
        touch ${prefix}.2.sai
        """
    }
}
